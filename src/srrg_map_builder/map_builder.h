#pragma once
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <string>
#include <queue>
#include <map>
#include <unordered_map>
#include <boost/bimap.hpp>
#include <Eigen/Core>
#include "srrg_core_map/map_node.h"
#include "srrg_core_map/map_node_list.h"
#include "srrg_core_map/local_map_with_traversability.h"


namespace srrg_map_builder{

typedef std::map<Eigen::Vector2f,srrg_core_map::MapNode*> Vector2fMapNodeMap;

class AxisAlignedSquare {

public:
    AxisAlignedSquare(){}

    AxisAlignedSquare(Eigen::Vector2f min, Eigen::Vector2f max):_min(min),_max(max){
        _halfWidth = (_max.x() - _min.x())/2;
        _center << _min.x() + _halfWidth,
                _min.y() + _halfWidth;
    }

    AxisAlignedSquare(Eigen::Vector2f center, float halfWidth):_center(center),_halfWidth(halfWidth){
        _min << _center.x() - _halfWidth,
                _center.y() - _halfWidth;
        _max << _center.x() + _halfWidth,
                _center.y() + _halfWidth;
    }

    inline void splitSquare(AxisAlignedSquare *squares){
        int idx = 0;
        float halfWidth = _halfWidth/2;
        for(int j = -1; j < 2; j += 2)
            for(int i = -1; i < 2; i += 2){
                squares[idx] = AxisAlignedSquare(Eigen::Vector2f (_center.x()+i*halfWidth,
                                                                  _center.y()+j*halfWidth),
                                                 halfWidth);
                idx++;
            }
    }

    inline bool inRange(const Eigen::Vector2f& point){
        return (point.x() >= _min.x() && point.x() <= _max.x() &&
                point.y() >= _min.y() && point.y() <= _max.y()) ? true : false;
    }

    inline const Eigen::Vector2f& min() const {return _min;}
    inline const Eigen::Vector2f& max() const {return _max;}
    inline const Eigen::Vector2f& center() const {return _center;}

protected:
    Eigen::Vector2f _min;
    Eigen::Vector2f _max;
    Eigen::Vector2f _center;
    float _halfWidth;
};

class QuadtreeNode {

public:
    QuadtreeNode(int depth, std::vector<Eigen::Vector2f> points, Eigen::Vector2f min, Eigen::Vector2f max):_depth(depth){
        if(_depth < 0)
            return;
        else{
            _square = AxisAlignedSquare(min,max);
            _points = points;
            AxisAlignedSquare squares[4];
            _square.splitSquare(squares);
            std::vector<Eigen::Vector2f> splitted_points[4];

            for(int ii=0; ii < _points.size(); ii++)
                for(int id = 0; id < 4; id++)
                    if(squares[id].inRange(_points.at(ii)))
                        splitted_points[id].push_back(_points.at(ii));

            for(int i = 0; i < 4; i++)
                if(!splitted_points[i].empty())
                    _children[i] = new QuadtreeNode(depth-1,splitted_points[i],squares[i].min(),squares[i].max());

        }
    }

    inline QuadtreeNode* getChild1(){return _children[0];}
    inline QuadtreeNode* getChild2(){return _children[1];}
    inline QuadtreeNode* getChild3(){return _children[2];}
    inline QuadtreeNode* getChild4(){return _children[3];}
    inline const int depth() const {return _depth;}
    inline const AxisAlignedSquare& square() const {return _square;}
    inline const std::vector<Eigen::Vector2f>& points() const {return _points;}

private:
    int _depth;
    AxisAlignedSquare _square;
    std::vector<Eigen::Vector2f> _points;
    QuadtreeNode* _children[4] = {NULL, NULL, NULL, NULL};
};


class Merger {

public:
    Merger(int depth_ = 0, float resolution_ = 0.05, float range_ = 1):
        _depth(depth_),
        _resolution(resolution_),
        _range(range_),
        _merged_local_maps(new srrg_core_map::MapNodeList),
        _quadtree(NULL) {
        _count = 0;
    }

    void computeBoundingBox(srrg_core_map::MapNodeList *local_maps_);

    inline void buildQuadtree(){
        std::cerr << "Building quadtree..." << std::endl;
        _quadtree = new QuadtreeNode(_depth,_points,_min,_max);
    }

    inline srrg_core_map::MapNodeList* execute(){
        std::cerr << "Executing merging algorithm..." << std::endl;
        visit(_quadtree);
        std::cerr << "Created " << _count << " local maps!!!" << std::endl;
        return _merged_local_maps;
    }

    inline void visualizeQuadtree(){
        _image = cv::Mat (_size.x(),_size.y(),CV_8UC3,cv::Scalar(255,255,255));
        visualize(_quadtree);
        cv::namedWindow("Image",cv::WINDOW_NORMAL);
        cv::imshow("Image",_image);
        cv::waitKey(0);
    }

private:
    int _count;
    int _depth;
    float _resolution;
    float _range;
    Vector2fMapNodeMap _local_maps;
    std::vector<Eigen::Vector2f> _points;
    srrg_core_map::MapNodeList* _merged_local_maps;
    Eigen::Vector2f _min;
    Eigen::Vector2f _max;
    Eigen::Vector2i _size;
    QuadtreeNode* _quadtree;
    cv::Mat _image;
    void visit(QuadtreeNode* quadtree);
    void visualize(QuadtreeNode *quadtree);
};

class Linker{

public:
    Linker(float distance_threshold_ = 5, float connectivity_threshold_ = 0.01, float resolution_ = 0.05):
        _distance_threshold(distance_threshold_),
        _connectivity_threshold(connectivity_threshold_),
        _resolution(resolution_),
        _local_maps(new srrg_core_map::MapNodeList),
        _relations(new srrg_core_map::BinaryNodeRelationSet){}

    inline void setInput(srrg_core_map::MapNodeList* local_maps_) {_local_maps = local_maps_;}

    srrg_core_map::BinaryNodeRelationSet* execute();

    bool addEdge(srrg_core_map::LocalMapWithTraversability* lmap1, srrg_core_map::LocalMapWithTraversability* lmap2);

    inline bool same(srrg_core_map::LocalMapWithTraversability* lmap1, srrg_core_map::LocalMapWithTraversability* lmap2){
        if(lmap1==lmap2)
            std::cerr << "Local maps are the same!" << std::endl;
        return (lmap1 == lmap2) ? true : false;
    }

    inline bool closeEnough(srrg_core_map::LocalMapWithTraversability* lmap1, srrg_core_map::LocalMapWithTraversability* lmap2){
        if((lmap1->transform().translation() - lmap2->transform().translation()).norm() > _distance_threshold)
            std::cerr << "Local maps are not close enough!" << std::endl;
        return ((lmap1->transform().translation() - lmap2->transform().translation()).norm() <= _distance_threshold) ? true : false;
    }

    inline bool alreadyConnected(srrg_core_map::LocalMapWithTraversability* lmap1, srrg_core_map::LocalMapWithTraversability* lmap2){
//        for(srrg_core_map::BinaryNodeRelationSet::iterator kt = _relations->begin(); kt != _relations->end(); kt++)
//            if((*kt)->from()->getId() == lmap1->getId() && (*kt)->to()->getId() == lmap2->getId() ||
//                    (*kt)->from()->getId() == lmap2->getId() && (*kt)->to()->getId() == lmap1->getId() ) {
//                std::cerr << "Local maps are already connected!" << std::endl;
//                return true;
//            }
//        return false;
        for(srrg_core_map::BinaryNodeRelationSet::iterator kt = _relations->begin(); kt != _relations->end(); kt++)
            if((*kt)->from() == lmap1 && (*kt)->to() == lmap2 ||
                    (*kt)->from() == lmap2 && (*kt)->to() == lmap1) {
                std::cerr << "Local maps are already connected!" << std::endl;
                return true;
            }
        return false;
    }

protected:
    float _distance_threshold;
    float _connectivity_threshold;
    float _resolution;
    srrg_core_map::MapNodeList* _local_maps;
    srrg_core_map::BinaryNodeRelationSet* _relations;
};
}
