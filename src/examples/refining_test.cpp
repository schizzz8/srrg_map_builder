#include <fstream>
#include <stdexcept>
#include <iostream>
#include <math.h>
#include <limits>
#include <queue>
#include <ctime>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <utility>
#include <map>
#include <functional>
#include <string>
#include <unordered_map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <srrg_boss/deserializer.h>
#include <srrg_boss/serializer.h>
#include <srrg_boss/trusted_loaders.h>
#include <srrg_system_utils/system_utils.h>
#include "srrg_core_map/image_map_node.h"
#include "srrg_core_map/local_map.h"
#include "srrg_core_map/local_map_with_traversability.h"

using namespace std;
using namespace srrg_boss;
using namespace srrg_core;
using namespace srrg_core_map;

namespace std {
template<>
bool std::less<Eigen::Vector3f>::operator ()(const Eigen::Vector3f& a,const Eigen::Vector3f& b) const {
    for(size_t i=0;i<3;++i) {
        if(a[i]<b[i]) return true;
        if(a[i]>b[i]) return false;
    }
    return false;
}
}

typedef std::map<Eigen::Vector3f,int> PointIndexMap;

class Cell {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Cell(const Eigen::Vector3i& idx=Eigen::Vector3i::Zero()):
        _idx(idx){
        _ground = false;
        _first_cloud = false;
        _second_cloud = false;
    }

    inline bool operator < (const Cell& c) const {
        for (int i=0; i<3; i++){
            if (_idx[i]<c._idx[i])
                return true;
            if (_idx[i]>c._idx[i])
                return false;
        }
        return false;
    }

    inline bool operator == (const Cell& c) const {
        for (int i=0; i<3; i++)
            if(_idx[i] != c._idx[i])
                return false;
        return true;
    }

    inline void setCenter(Eigen::Vector3f origin, float resolution) {
        _center = origin + _idx.cast<float>()*resolution + Eigen::Vector3f(resolution/2,resolution/2,resolution/2);
    }

    inline const bool ground() const {return _ground;}
    inline void setGround(const bool& ground_){_ground = ground_;}
    inline const bool overlap() const {return _first_cloud & _second_cloud;}
    //inline void setOverlap(const bool& overlap_){_overlap = overlap_;}

    Eigen::Vector3i _idx;
    Eigen::Vector3f _center;
    PointIndexMap _points;
    bool _ground;
    bool _first_cloud,_second_cloud;

};

template<typename T>
struct matrix_hash : std::unary_function<T, size_t> {
    std::size_t operator()(T const& matrix) const {
        size_t seed = 0;
        for (size_t i = 0; i < matrix.size(); ++i) {
            auto elem = *(matrix.data() + i);
            seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
typedef std::unordered_map<Eigen::Vector3i,Cell*,matrix_hash<Eigen::Vector3i> > Vector3iCellPtrMap;

class SparseGrid : public Vector3iCellPtrMap {
public:
    SparseGrid (float resolution_ = 0.05,
                Eigen::Vector3f origin_ = Eigen::Vector3f::Zero(),
                Eigen::Vector3i dimensions_ = Eigen::Vector3i::Zero(),
                int half_ = std::numeric_limits<int>::max(),
                float threshold_ = 0.01)
        :_resolution(resolution_),
          _origin(origin_),
          _dimensions(dimensions_),
          _half(half_),
          _threshold(threshold_){
        _inverse_resolution = 1./_resolution;
        _num_cells = _dimensions.x()*_dimensions.y()*_dimensions.z();
        _cloud = new Cloud;
        _traversability_map = new TraversabilityMap;
    }

    void insertCloud(const srrg_core_map::Cloud& cloud){
        for(size_t ii=0; ii < cloud.size(); ii++) {

            Eigen::Vector3i idx = toGrid(cloud.at(ii).point());

            if(hasCell(idx)) {
                at(idx)->_points.insert(std::pair<Eigen::Vector3f,int> (cloud.at(ii).point(),ii));

                if(ii < _half)
                    at(idx)->_first_cloud = true;
                else
                    at(idx)->_second_cloud = true;
                //                if(at(idx)->_overlap == false){
                //                    bool first_cloud = false, second_cloud = false;
                //                    for(PointIndexMap::iterator it = at(idx)->_points.begin();
                //                        it != at(idx)->_points.end();
                //                        it++){
                //                        if(it->second < _half)
                //                            first_cloud = true;
                //                        else
                //                            second_cloud = true;
                //                    }
                //                    if(first_cloud && second_cloud)
                //                        at(idx)->setOverlap(true);
                //                }

            }
            else {
                Cell* cell = new Cell(idx);
                cell->setCenter(_origin,_resolution);
                cell->_points.insert(std::pair<Eigen::Vector3f,int> (cloud.at(ii).point(),ii));
                Vector3iCellPtrMap::iterator it = begin();
                insert(it,std::pair<Eigen::Vector3i,Cell*>(idx,cell));
            }
        }

    }

    srrg_core_map::Cloud* extractCloud(){

        for(Vector3iCellPtrMap::iterator it = begin();
            it != end();
            it++){

            const PointIndexMap& points = it->second->_points;
            Eigen::Vector3f centroid = Eigen::Vector3f::Zero();

            for(PointIndexMap::const_iterator jt = points.begin();
                jt != points.end();
                jt++){
                const Eigen::Vector3f& point = jt->first;
                centroid += point;
            }

            centroid /= points.size();
            _cloud->push_back(RichPoint(centroid));
        }

        return _cloud;
    }

    srrg_core_map::TraversabilityMap* extractSurface(){
        IntImage indices;
        FloatImage elevations;
        UnsignedCharImage classified;
        Eigen::Vector3f bottom;

        indices.release();
        elevations.release();

        bottom = _origin;
        bottom.z() = 0;

        int cols = _dimensions.x();
        int rows = _dimensions.y();

        indices.create(rows,cols);
        elevations.create(rows,cols);
        classified.create(rows,cols);

        indices = -1;
        elevations = 5;
        classified = 127;

        float robot_climb_step=0.1;
        float robot_height=0.5;

        for (int i = 0; i < _cloud->size(); i++){
            const RichPoint& point = _cloud->at(i);
            float z = point.point().z();
            Eigen::Vector3f projected_point = (point.point() - bottom)*_inverse_resolution;
            int row = projected_point.y();
            int col = projected_point.x();
            if(row>=indices.rows || row<0)
                continue;
            if(col>=indices.cols || col<0)
                continue;
            float& height = elevations.at<float> (row,col);
            int& idx = indices.at<int> (row,col);
            if(z<height){
                height = z;
                idx = i;
            }
        }

        for (int i = 0; i < _cloud->size(); i++){
            const RichPoint& point = _cloud->at(i);
            float z = point.point().z();
            Eigen::Vector3f projected_point = (point.point() - bottom)*_inverse_resolution;
            int row = projected_point.y();
            int col = projected_point.x();
            if(row>=indices.rows || row<0)
                continue;
            if(col>=indices.cols || col<0)
                continue;
            float& height = elevations.at<float> (row,col);
            int& idx = indices.at<int> (row,col);
            float min_obstacle_height = height+robot_climb_step;
            float max_obstacle_height = height+robot_height;
            if (z < min_obstacle_height)
                continue;
            if (z > max_obstacle_height)
                continue;
            idx = -2;
            height = z;
        }

        for (int r=0; r<indices.rows; r++)
            for (int c=0; c<indices.cols; c++) {
                int idx = indices.at<int>(r,c);
                if (idx==-1)
                    continue;
                if (idx<-1){
                    classified.at<unsigned char>(r,c)=255;
                    continue;
                }
                classified.at<unsigned char>(r,c)=0;
            }

        for (int r=1; r<classified.rows-1; r++)
            for (int c=1; c<classified.cols-1; c++) {
                unsigned char & cell=classified.at<unsigned char>(r,c);
                if (cell!=255)
                    continue;
                bool one_big=false;
                for (int rr=-1; rr<=1; rr++)
                    for (int cc=-1; cc<=1; cc++) {
                        if (rr==0 && cc==0)
                            continue;
                        one_big |= classified.at<unsigned char>(r+rr,c+cc)==255;
                    }
                if (! one_big) {
                    cell=0;
                }
            }

        //        cv::imshow("map",classified);
        //        cv::waitKey(0);

        _traversability_map = new TraversabilityMap (classified);

        //        for (int i = 0; i < _cloud->size(); i++){
        //            const RichPoint& point = _cloud->at(i);
        //            Eigen::Vector3i idx = toGrid(point.point());
        //            unsigned char & cell=classified.at<unsigned char>(idx.y(),idx.x());
        //            if(cell == 0)
        //                at(idx)->setGround(true);
        //        }
        for (int r=0; r<indices.rows; r++)
            for (int c=0; c<indices.cols; c++) {
                int idx = indices.at<int>(r,c);
                if (idx<0)
                    continue;
                unsigned char & cell=classified.at<unsigned char>(r,c);
                if(cell==0){
                    const RichPoint& point = _cloud->at(idx);
                    Eigen::Vector3i idx = toGrid(point.point());
                    at(idx)->setGround(true);
                }
            }


        return _traversability_map;

    }

    bool checkConnectivity(){
        float ground_count=0;
        float overlap_count=0;
        for(Vector3iCellPtrMap::iterator it = begin();
            it != end();
            it++){
            Cell* cell = it->second;
            if(cell->ground() == true){
                ground_count++;
                if(cell->overlap() == true)
                    overlap_count++;
            }
        }
        //        cerr << "Overlap: " << overlap_count << endl;
        //        cerr << "Ground: " << ground_count << endl;
        //        cerr << "Connection percentage: " << overlap_count/ground_count << endl;
        return (overlap_count/ground_count > _threshold) ? true : false;
    }

    inline float resolution(){ return _resolution;}
    inline const Eigen::Vector3i dimensions(){ return _dimensions;}
    inline const Eigen::Vector3f origin(){ return _origin;}
    inline int numCells(){ return _num_cells;}

    inline const Eigen::Vector3i toGrid(const Eigen::Vector3f& point) const {
        return ((point - _origin)*_inverse_resolution).cast<int>();
    }
    inline const Eigen::Vector3f toWorld(const Eigen::Vector3i& cell) const{
        return (_origin + cell.cast<float>() *_resolution);
    }


protected:
    float _resolution;
    float _inverse_resolution;
    Eigen::Vector3f _origin;
    Eigen::Vector3i _dimensions;
    int _num_cells;
    int _half;
    float _threshold;
    Cloud* _cloud;
    TraversabilityMap* _traversability_map;

    inline const bool hasCell(const Eigen::Vector3i& idx){
        Vector3iCellPtrMap::iterator it = find(idx);
        return (it != end()) ? true : false;
    }
};

bool same(LocalMap* lmap1, LocalMap* lmap2){
    return (lmap1->getId() == lmap2->getId()) ? true : false;
}

bool closeEnough(LocalMap* lmap1, LocalMap* lmap2, float distance_threshold){
    return ((lmap1->transform().translation() - lmap2->transform().translation()).norm() <= distance_threshold) ? true : false;
}

bool alreadyConnected(const BinaryNodeRelationSet* relations,LocalMap* lmap1, LocalMap* lmap2){
    for(BinaryNodeRelationSet::iterator kt = relations->begin(); kt != relations->end(); kt++)
        if((*kt)->from()->getId() == lmap1->getId() && (*kt)->to()->getId() == lmap2->getId() ||
                (*kt)->from()->getId() == lmap2->getId() && (*kt)->to()->getId() == lmap1->getId() ) {
            return true;
        }
    return false;
}

bool addEdge(LocalMap* lmap1, LocalMap* lmap2){

    float resolution = 0.01;
    Eigen::Vector3f origin = Eigen::Vector3f::Zero();
    Eigen::Vector3i dimensions = Eigen::Vector3i::Zero();

    Cloud* cloud = new Cloud;
    Cloud transformed_cloud;
    lmap1->cloud()->transform(transformed_cloud,lmap1->transform());
    cloud->add(transformed_cloud);
    int size1 = cloud->size();
    transformed_cloud.clear();
    lmap2->cloud()->transform(transformed_cloud,lmap2->transform());
    cloud->add(transformed_cloud);

    Eigen::Vector3f lower,higher;
    cloud->computeBoundingBox(lower,higher);
    origin = lower - Eigen::Vector3f(5*resolution,5*resolution,5*resolution);
    dimensions = (((higher + Eigen::Vector3f(5*resolution,5*resolution,5*resolution))-lower)/resolution).cast<int> ();
    //    cerr << "Grid size: " << dimensions.transpose() << endl;
    //    cerr << "Half: " << size1 << endl;
    SparseGrid grid (resolution,origin,dimensions,size1);
    grid.insertCloud(*cloud);
    grid.extractCloud();
    grid.extractSurface();

    return grid.checkConnectivity();

}



// Help objects to force linking
BaseCameraInfo cinfo;
ImageMapNode tnode;
LocalMap lmap;
BinaryNodeRelation rel;



const char* banner[] = {
    "srrg_trajectory_loader_app: example on how to load a set of boss objects constituting a boss map",
    "usage:",
    " srrg_trajectory_loader_app <boss log file>",
    0
};

int main(int argc, char** argv) {
    if(argc < 2 || !strcmp(argv[1],"-h")) {
        printBanner(banner);
        return 0 ;
    }

    string input_filename = argv[1];

    std::list<Serializable*> objects;
    MapNodeList* lmaps = new MapNodeList;
    BinaryNodeRelationSet* relations = new BinaryNodeRelationSet;
    Deserializer des;
    des.setFilePath(input_filename);
    Serializable* o;
    while((o = des.readObject())) {
        LocalMap* lmap = dynamic_cast<LocalMap*>(o);
        if (lmap){
            cerr << "Local map: " << lmap->getId() << endl;
            cerr << "\t>> " << lmap->cloud()->size() << " points" << endl;
            cerr << "\t>> " << lmap->nodes().size() << " nodes" << endl;
            cerr << "\t>> " << lmap->relations().size() << " relations" << endl;
            lmaps->addElement(lmap);
        }
        objects.push_back(o);
    }

    float distance_threshold = 5;
    int count = 0;

    for(MapNodeList::iterator it = lmaps->begin(); it != lmaps->end(); it++) {
        LocalMap* lmap1 = dynamic_cast<LocalMap*> (it->get());
        for(MapNodeList::iterator jt = lmaps->begin(); jt != lmaps->end(); jt++){
            LocalMap* lmap2 = dynamic_cast<LocalMap*> (jt->get());
            //cerr << "Checking connectivity for lmap " << lmap1->getId() << " and lmap " << lmap2->getId() << endl;
            if(!same(lmap1,lmap2) && closeEnough(lmap1,lmap2,distance_threshold) && !alreadyConnected(relations,lmap1,lmap2))
                if(addEdge(lmap1,lmap2)) {
                    BinaryNodeRelation* rel = new BinaryNodeRelation(lmap1,lmap2,lmap1->transform().inverse()*lmap2->transform());
                    relations->insert(std::tr1::shared_ptr<BinaryNodeRelation>(rel));
                    count++;
                    //cerr << "adding edge!!!" << endl << endl;
                }
            //                    else
            //                    cerr << "local maps are not connected!" << endl << endl;
        }
    }
    cerr << "Added " << count << " edges" << endl;

    for(BinaryNodeRelationSet::iterator it = relations->begin(); it != relations->end(); it++) {
        BinaryNodeRelation* relation = dynamic_cast<BinaryNodeRelation*> (it->get());
        if(relation)
            objects.push_back(relation);
    }

    cerr << "Writing: " << objects.size() << " objects" << endl;
    Serializer ser;
    string output_filename = input_filename.substr(0,input_filename.find("."))+"_connected.lmaps";
    ser.setFilePath(output_filename);
    ser.setBinaryPath(output_filename + ".d/<classname>.<nameAttribute>.<id>.<ext>");
    for(std::list<Serializable*>::iterator it = objects.begin(); it != objects.end(); it++){
        Serializable* s = *it;
        ser.writeObject(*s);
    }

    return 0;
}
