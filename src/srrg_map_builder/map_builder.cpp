#include "map_builder.h"


using namespace std;

namespace std {
template<>
bool std::less<Eigen::Vector2f>::operator ()(const Eigen::Vector2f& a,const Eigen::Vector2f& b) const {
    for(size_t i=0;i<2;++i) {
        if(a[i]<b[i]) return true;
        if(a[i]>b[i]) return false;
    }
    return false;
}
}

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

namespace srrg_map_builder {

using namespace srrg_core;
using namespace srrg_core_map;

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
                int half_ = std::numeric_limits<int>::max())
        :_resolution(resolution_),
          _origin(origin_),
          _dimensions(dimensions_),
          _half(half_){
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
            } else {
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

        _traversability_map = new TraversabilityMap (classified);

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

    bool checkConnectivity(float connectivity_threshold){
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
        cerr << "Overlap: " << overlap_count << endl;
        cerr << "Ground: " << ground_count << endl;
        cerr << "Connection percentage: " << overlap_count/ground_count << endl;
        return (overlap_count/ground_count > connectivity_threshold) ? true : false;
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
    Cloud* _cloud;
    TraversabilityMap* _traversability_map;

    inline const bool hasCell(const Eigen::Vector3i& idx){
        Vector3iCellPtrMap::iterator it = find(idx);
        return (it != end()) ? true : false;
    }
};

void Merger::computeBoundingBox(MapNodeList* local_maps_){
    int n = 0;
    Eigen::Vector2f centroid = Eigen::Vector2f::Zero();
    Eigen::Vector2f max,min;
    float xmin=std::numeric_limits<float>::max();
    float xmax=std::numeric_limits<float>::min();
    float ymin=std::numeric_limits<float>::max();
    float ymax=std::numeric_limits<float>::min();
    for(MapNodeList::iterator it = local_maps_->begin(); it != local_maps_->end(); it++){
        LocalMap* lmap = dynamic_cast<LocalMap*> (it->get());
        Eigen::Vector2f point = Eigen::Vector2f (lmap->transform().translation().x(),
                                                 lmap->transform().translation().y());
        _points.push_back(point);
        _local_maps.insert(std::pair<Eigen::Vector2f,MapNode*>(point,it->get()));
        xmax = xmax > point.x()+_range ? xmax : point.x()+_range;
        ymax = ymax > point.y()+_range ? ymax : point.y()+_range;
        xmin = xmin < point.x()-_range ? xmin : point.x()-_range;
        ymin = ymin < point.y()-_range ? ymin : point.y()-_range;
        centroid += point;
        n++;
    }
    centroid /= n;
    min = Eigen::Vector2f (xmin,ymin);
    max = Eigen::Vector2f (xmax,ymax);

    float radius = (max-centroid).norm() + 1*_resolution;
    _min = Eigen::Vector2f(centroid.x()-radius,centroid.y()-radius);
    _max = Eigen::Vector2f(centroid.x()+radius,centroid.y()+radius);

    cerr << endl << "Bounding Box: " << endl;
    cerr << "\t>>Lower: " << _min.transpose() << endl;
    cerr << "\t>>Upper: " << _max.transpose() << endl;

    _size = Eigen::Vector2i ((_max.x()-_min.x())/_resolution,
                             (_max.x()-_min.x())/_resolution);

    cerr << "Grid size: " << _size.x() << "x" << _size.y() << endl;

}

void Merger::visualize(QuadtreeNode *quadtree){
    if(quadtree){
        if(quadtree->depth() == 0){

            float xmin = quadtree->square().min().x();
            float ymin = quadtree->square().min().y();
            float xmax = quadtree->square().max().x();
            float ymax = quadtree->square().max().y();
            cv::rectangle(_image,
                          cv::Point ((xmin-_min.x())/_resolution,(ymin-_min.y())/_resolution),
                          cv::Point ((xmax-_min.x())/_resolution,(ymax-_min.y())/_resolution),
                          cv::Scalar (0,0,255));

            for(vector<Eigen::Vector2f>::const_iterator it = quadtree->points().begin();
                it != quadtree->points().end();
                it++){
                Eigen::Vector2f point = *it;
                cv::circle(_image,
                           cv::Point ((point.x()-_min.x())/_resolution,(point.y()-_min.y())/_resolution),
                           1,
                           cv::Scalar (255,0,0));
            }

        }
        visualize(quadtree->getChild1());
        visualize(quadtree->getChild2());
        visualize(quadtree->getChild3());
        visualize(quadtree->getChild4());
    } else {
        return;
    }
}

void Merger::visit(QuadtreeNode* quadtree){
    if(quadtree){
        if(quadtree->depth() == 0) {
            _count++;

            MapNodeList* lmaps = new MapNodeList;
            for(vector<Eigen::Vector2f>::const_iterator it = quadtree->points().begin();
                it != quadtree->points().end();
                it++){
                Vector2fMapNodeMap::iterator iter = _local_maps.find(*it);
                if(iter != _local_maps.end())
                    lmaps->addElement(iter->second);
            }

            Eigen::Isometry3f transform = Eigen::Isometry3f::Identity();
            transform = lmaps->middlePose();
            LocalMapWithTraversability* reference = new LocalMapWithTraversability (transform);
            reference->setCloud(new Cloud());
            Cloud transformed_cloud;
            Eigen::Vector3f lower,higher;
            Eigen::Vector3f origin = Eigen::Vector3f::Zero();
            Eigen::Vector3i dimensions = Eigen::Vector3i::Zero();

            for(MapNodeList::iterator it = lmaps->begin(); it != lmaps->end(); it++){
                LocalMap* current = dynamic_cast<LocalMap*> (it->get());
                if(current){
                    Cloud().swap(transformed_cloud);
                    current->cloud()->transform(transformed_cloud,transform.inverse()*current->transform());
                    transformed_cloud.add(*(reference->cloud()));
                    transformed_cloud.computeBoundingBox(lower,higher);
                    origin = lower - Eigen::Vector3f(5*_resolution,5*_resolution,5*_resolution);
                    dimensions = (((higher + Eigen::Vector3f(5*_resolution,5*_resolution,5*_resolution))-lower)/_resolution).cast<int> ();
                    SparseGrid grid (_resolution,origin,dimensions);
                    grid.insertCloud(transformed_cloud);
                    reference->setCloud(grid.extractCloud());
                    reference->setTraversabilityMap(grid.extractSurface());
                    reference->setLower(origin);
                    reference->setUpper(higher);
                    reference->setResolution(_resolution);
                }
            }
            _merged_local_maps->addElement(reference);
        }
        visit(quadtree->getChild1());
        visit(quadtree->getChild2());
        visit(quadtree->getChild3());
        visit(quadtree->getChild4());
    }
    else
        return;
}

bool Linker::addEdge(LocalMapWithTraversability* lmap1, LocalMapWithTraversability* lmap2){
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
    origin = lower - Eigen::Vector3f(5*_resolution,5*_resolution,5*_resolution);
    dimensions = (((higher + Eigen::Vector3f(5*_resolution,5*_resolution,5*_resolution))-lower)/_resolution).cast<int> ();
    SparseGrid grid (_resolution,origin,dimensions,size1);
    grid.insertCloud(*cloud);
    grid.extractCloud();
    grid.extractSurface();

    return grid.checkConnectivity(_connectivity_threshold);

}

BinaryNodeRelationSet *Linker::execute(){
    for(MapNodeList::iterator it = _local_maps->begin(); it != _local_maps->end(); it++) {
        LocalMapWithTraversability* lmap1 = dynamic_cast<LocalMapWithTraversability*> (it->get());
        for(MapNodeList::iterator jt = _local_maps->begin(); jt != _local_maps->end(); jt++){
            LocalMapWithTraversability* lmap2 = dynamic_cast<LocalMapWithTraversability*> (jt->get());
            cerr << endl << "Checking connectivity for lmap " << lmap1 << " and lmap " << lmap2 << endl;
            if(!same(lmap1,lmap2) && closeEnough(lmap1,lmap2) && !alreadyConnected(lmap1,lmap2))
                if(addEdge(lmap1,lmap2)) {
                    BinaryNodeRelation* rel = new BinaryNodeRelation(lmap1,lmap2,lmap1->transform().inverse()*lmap2->transform());
                    _relations->insert(std::tr1::shared_ptr<BinaryNodeRelation>(rel));
                }
            cerr << "#######################################################################################" << endl;
        }
    }
    cerr << "Added " << _relations->size() << " edges" << endl;
    return _relations;
}

}
