/**
 * @file softrenderer.cpp
 * @author Sk. Mohammadul Haque
 * @version 0.1.0.0
 * @copyright
 * Copyright (c) 2019 Sk. Mohammadul Haque.
 * @brief This header file contains definitions of all functions of softrenderer.
 */
#define SOFTRENDERER_API_EXPORTS
#include "../include/softrenderer.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <stdint.h>
#include <algorithm>
#include <assert.h>

/** \cond HIDDEN_SYMBOLS */
#undef max
using namespace std;
template <typename T> int sgn(T val)
{
    return (T(0)<val)-(val<T(0));
}

#ifndef isnan
#define isnan(x) ((x)!=(x))
#endif

template <typename T>
class Vertex
{
public:
    T x, y;
    Vertex(T x_, T y_):x(x_), y(y_) {}
    Vertex():x(0.0), y(0.0) {}
};

template <typename T>
using VVert = std::vector<Vertex<T>>;

template <typename T>
class Triangle
{
public:
    Vertex<T> vt1, vt2, vt3;
    VVert<T> vv;
    bool is_culled;
    Triangle(T x1, T y1, T x2, T y2, T x3, T y3):vt1(x1,y1), vt2(x2,y2), vt3(x3,y3), is_culled(false) {}
    void sort_vertices_ascending_by_y();
    void fill_bottom_flat_triangle(Vertex<T> &v1, Vertex<T> &v2, Vertex<T> &v3);
    void fill_top_flat_triangle(Vertex<T> &v1, Vertex<T> &v2, Vertex<T> &v3);
    void draw_triangle();
    void check();
};

template <typename T>
void Triangle<T>::sort_vertices_ascending_by_y()
{
    Vertex<T> vTmp;
    if(vt1.y>vt2.y)
    {
        vTmp = vt1;
        vt1 = vt2;
        vt2 = vTmp;
    }
    if(vt1.y>vt3.y)
    {
        vTmp = vt1;
        vt1 = vt3;
        vt3 = vTmp;
    }
    if(vt2.y>vt3.y)
    {
        vTmp = vt2;
        vt2 = vt3;
        vt3 = vTmp;
    }
}

template <typename T>
void Triangle<T>::check()
{
    Vertex<T> v21(vt2.x-vt1.x, vt2.y-vt1.y);
    Vertex<T> v32(vt3.x-vt2.x, vt3.y-vt2.y);
    is_culled = v21.x*v32.y-v21.y*v32.x>0;
}

template <typename T>
void Triangle<T>::fill_bottom_flat_triangle(Vertex<T> &v1, Vertex<T> &v2, Vertex<T> &v3)
{
    T invslope1;
    T invslope2;
    T curx1 = v1.x;
    T curx2 = v1.x;
    T mn, mx;
    if(v3.x>v2.x)
    {
        invslope1 = (v2.x-v1.x)/(v2.y-v1.y);
        invslope2 = (v3.x-v1.x)/(v3.y-v1.y);
        mn = std::min(v2.x, v1.x);
        mx = std::max(v3.x, v1.x);
    }
    else
    {
        invslope2 = (v2.x-v1.x)/(v2.y-v1.y);
        invslope1 = (v3.x-v1.x)/(v3.y-v1.y);
        mn = std::min(v3.x, v1.x);
        mx = std::max(v2.x, v1.x);
    }
#ifdef __MY_DEBUG__
    CHECK_LINE_Y(std::floor(v1.y), std::ceil(v2.y), <)
#endif
    for(T scanlineY=std::floor(v1.y); scanlineY<=std::ceil(v2.y); ++scanlineY)
    {
#ifdef __MY_DEBUG__
        CHECK_LINE_X(floor(v1.y), ceil(v2.y), <)
#endif
        for(T scanlineX=std::floor(std::max(curx1, mn)); scanlineX<=std::ceil(std::min(curx2, mx)); ++scanlineX)
        {
            vv.push_back(Vertex<T>(scanlineX,scanlineY));
        }
        curx1 += invslope1;
        curx2 += invslope2;
    }
}

template <typename T>
void Triangle<T>::fill_top_flat_triangle(Vertex<T> &v1, Vertex<T> &v2, Vertex<T> &v3)
{
    T invslope1;
    T invslope2;
    T curx1 = v3.x;
    T curx2 = v3.x;
    T mn, mx;
    if(v2.x>v1.x)
    {
        invslope1 = (v3.x-v1.x)/(v3.y-v1.y);
        invslope2 = (v3.x-v2.x)/(v3.y-v2.y);
        mn = std::min(v1.x, v3.x);
        mx = std::max(v2.x, v3.x);
    }
    else
    {
        invslope2 = (v3.x-v1.x)/(v3.y-v1.y);
        invslope1 = (v3.x-v2.x)/(v3.y-v2.y);
        mn = std::min(v2.x, v3.x);
        mx = std::max(v1.x, v3.x);
    }
    for(T scanlineY=std::ceil(v3.y); scanlineY>=std::floor(v1.y); --scanlineY)
    {
        curx1 -= invslope1;
        curx2 -= invslope2;
        for(T scanlineX=std::floor(std::max(curx1, mn)); scanlineX<=std::ceil(std::min(curx2, mx)); ++scanlineX)
        {
            vv.push_back(Vertex<T>(scanlineX,scanlineY));
        }
    }
}

template <typename T>
void Triangle<T>::draw_triangle()
{
    sort_vertices_ascending_by_y();
    if(vt2.y==vt3.y)
    {
        fill_bottom_flat_triangle(vt1, vt2, vt3);
    }
    else if(vt1.y==vt2.y)
    {
        fill_top_flat_triangle(vt1, vt2, vt3);
    }
    else
    {
        Vertex<T> v4((vt1.x+((vt2.y-vt1.y)/(vt3.y-vt1.y))*(vt3.x-vt1.x)), vt2.y);
        fill_bottom_flat_triangle(vt1, vt2, v4);
        fill_top_flat_triangle(vt2, v4, vt3);
    }
}

/** \endcond */

/** \brief Perspective-renders a given mesh to a depthmap and optionally scalarmap
 *
 * \param[in] m Input mesh (must contain faces, face-scalars are preferred)
 * \param[out] depthmap Column-major depth map
 * \param[out] scalarmap Column-major scalar map
 * \param[in] f Camera focal-length (in px)
 * \param[in] h Map height (in px)
 * \param[in] w Map width (in px)
 * \param[in] doscalar Render scalar map (true/false)
 * \param[in] cutoff Render front z-cut-off (default - 0.0)
 * \param[in] eps Render tolerance (default - 0.0001)
 * \return Error code (0-success)
 *
 */

SOFTRENDERER_API int sr_mesh_render_persp(MESH m, std::vector<FLOATDATA>& depthmap, std::vector<FLOATDATA>& scalarmap, FLOATDATA f, INTDATA h, INTDATA w, bool doscalar, FLOATDATA cutoff, FLOATDATA eps)
{
    INTDATA nv, nf, xx, yy;
    std::vector<FLOATDATA> x, y;
    FLOATDATA Tdist, fw, fh, a, b, c, d, px, py, pz; //, mn, mx;

    fw = w/2.0;
    fh = h/2.0;

    MESH_VERTEX mv = m->vertices;
    nv = m->num_vertices;
    MESH_SCALAR mvs = m->vscalars;
    MESH_FACE mf = m->faces;
    nf = m->num_faces;
    MESH_SCALAR mfs = m->fscalars;

    x.resize(nv);
    y.resize(nv);

    for(INTDATA i=0; i<nv; ++i)
    {
        const MESH_VERTEX mvc = mv+i;
        if(mvc->z<-cutoff)
        {
            x[i] = fw-f*mvc->x/mvc->z;
            y[i] = fh+f*mvc->y/mvc->z;
            //if(x[i]<0.0||y[i]<0.0||x[i]>=w||y[i]>=h) // this is valid only when sufficiently sampled
            if(x[i]<-w/4||y[i]<-h/4||x[i]>=(5*w)/4||y[i]>=(5*h)/4)
            {
                x[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
                y[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
            }
        }
        else
        {
            x[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
            y[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
        }
    }

    if(m->is_fscalars&&doscalar)
    {
        depthmap.resize(h*w);
        std::fill(depthmap.begin(), depthmap.end(), std::numeric_limits<FLOATDATA>::max());
        scalarmap.resize(h*w);
        std::fill(scalarmap.begin(), scalarmap.end(), std::numeric_limits<FLOATDATA>::quiet_NaN());

        for(INTDATA i=0; i<nf; ++i)
        {
            const INTDATA* mfcv = mf[i].vertices;
            const INTDATA vx1 = mfcv[0], vx2 = mfcv[1], vx3 = mfcv[2];
            if(isnan(x[vx1])||isnan(x[vx2])||isnan(x[vx3])) continue;
            Triangle<FLOATDATA> ct(x[vx1],y[vx1],x[vx2],y[vx2],x[vx3],y[vx3]);
            ct.check();
            ct.draw_triangle();
            a = (mv[vx2].y-mv[vx1].y)*(mv[vx3].z-mv[vx1].z)
                -(mv[vx3].y-mv[vx1].y)*(mv[vx2].z-mv[vx1].z);
            b = (mv[vx2].z-mv[vx1].z)*(mv[vx3].x-mv[vx1].x)
                -(mv[vx3].z-mv[vx1].z)*(mv[vx2].x-mv[vx1].x);
            c = (mv[vx2].x-mv[vx1].x)*(mv[vx3].y-mv[vx1].y)
                -(mv[vx3].x-mv[vx1].x)*(mv[vx2].y-mv[vx1].y);
            d = -(a*mv[vx1].x+b*mv[vx1].y+c*mv[vx1].z);
            px = -a/f;
            py = b/f;
            pz = a*fw/f-b*fh/f+c;
            FLOATDATA mn = -std::max(std::max(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            FLOATDATA mx = -std::min(std::min(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            for(size_t j=0; j<ct.vv.size(); ++j)
            {
                Tdist = d/(px*ct.vv[j].x+py*ct.vv[j].y+pz); // (made positive for in front of camera)
                if(Tdist<mn*0.995) Tdist = mn; // check here this is buggy
                if(Tdist>mx*1.005) Tdist = mx; // check here this is buggy
                if(Tdist<=0.0) continue;
                xx = static_cast<int>(ct.vv[j].x)-1;
                yy = static_cast<int>(ct.vv[j].y)-1;
                if(xx>=0 && xx<w && yy>=0 && yy<h)
                {
                    if (depthmap[xx*h + yy] > Tdist)
                    {
                        depthmap[xx*h + yy] = Tdist;
                        scalarmap[xx*h + yy] = mfs[i]; // assumed faces have scalars
                    }
                }
            }
        }
    }
    else if(m->is_vscalars&&doscalar)
    {
        depthmap.resize(h*w);
        std::fill(depthmap.begin(), depthmap.end(), std::numeric_limits<FLOATDATA>::max());
        scalarmap.resize(h*w);
        std::fill(scalarmap.begin(), scalarmap.end(), std::numeric_limits<FLOATDATA>::quiet_NaN());

        for(INTDATA i=0; i<nf; ++i)
        {
            const INTDATA* mfcv = mf[i].vertices;
            const INTDATA vx1 = mfcv[0], vx2 = mfcv[1], vx3 = mfcv[2];
            if(isnan(x[vx1])||isnan(x[vx2])||isnan(x[vx3])) continue;
            Triangle<FLOATDATA> ct(x[vx1],y[vx1],x[vx2],y[vx2],x[vx3],y[vx3]);
            ct.check();
            ct.draw_triangle();
            a = (mv[vx2].y-mv[vx1].y)*(mv[vx3].z-mv[vx1].z)
                -(mv[vx3].y-mv[vx1].y)*(mv[vx2].z-mv[vx1].z);
            b = (mv[vx2].z-mv[vx1].z)*(mv[vx3].x-mv[vx1].x)
                -(mv[vx3].z-mv[vx1].z)*(mv[vx2].x-mv[vx1].x);
            c = (mv[vx2].x-mv[vx1].x)*(mv[vx3].y-mv[vx1].y)
                -(mv[vx3].x-mv[vx1].x)*(mv[vx2].y-mv[vx1].y);
            d = -(a*mv[vx1].x+b*mv[vx1].y+c*mv[vx1].z);
            px = -a/f;
            py = b/f;
            pz = a*fw/f-b*fh/f+c;
            FLOATDATA mn = -std::max(std::max(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            FLOATDATA mx = -std::min(std::min(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            for(size_t j=0; j<ct.vv.size(); ++j)
            {
                Tdist = d/(px*ct.vv[j].x+py*ct.vv[j].y+pz); // (made positive for in front of camera)
                if(Tdist<mn*0.995) Tdist = mn; // check here this is buggy
                if(Tdist>mx*1.005) Tdist = mx; // check here this is buggy
                if(Tdist<=0.0) continue;
                xx = static_cast<int>(ct.vv[j].x)-1;
                yy = static_cast<int>(ct.vv[j].y)-1;
                if(xx>=0 && xx<w && yy>=0 && yy<h)
                {
                    if (depthmap[xx*h + yy] > Tdist)
                    {
                        depthmap[xx*h + yy] = Tdist;
                        scalarmap[xx*h + yy] = mvs[vx1];// assumed all vertex has same scalar
                    }
                }
            }
        }
    }
    else
    {
        depthmap.resize(h*w);
        std::fill(depthmap.begin(), depthmap.end(), std::numeric_limits<FLOATDATA>::max());

        for(INTDATA i=0; i<nf; ++i)
        {
            const INTDATA* mfcv = mf[i].vertices;
            const INTDATA vx1 = mfcv[0], vx2 = mfcv[1], vx3 = mfcv[2];
            if(isnan(x[vx1])||isnan(x[vx2])||isnan(x[vx3])) continue;
            Triangle<FLOATDATA> ct(x[vx1],y[vx1],x[vx2],y[vx2],x[vx3],y[vx3]);
            ct.check();
            ct.draw_triangle();
            a = (mv[vx2].y-mv[vx1].y)*(mv[vx3].z-mv[vx1].z)
                -(mv[vx3].y-mv[vx1].y)*(mv[vx2].z-mv[vx1].z);
            b = (mv[vx2].z-mv[vx1].z)*(mv[vx3].x-mv[vx1].x)
                -(mv[vx3].z-mv[vx1].z)*(mv[vx2].x-mv[vx1].x);
            c = (mv[vx2].x-mv[vx1].x)*(mv[vx3].y-mv[vx1].y)
                -(mv[vx3].x-mv[vx1].x)*(mv[vx2].y-mv[vx1].y);
            d = -(a*mv[vx1].x+b*mv[vx1].y+c*mv[vx1].z);
            px = -a/f;
            py = b/f;
            pz = a*fw/f-b*fh/f+c;
            FLOATDATA mn = -std::max(std::max(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            FLOATDATA mx = -std::min(std::min(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            for(size_t j=0; j<ct.vv.size(); ++j)
            {
                Tdist = d/(px*ct.vv[j].x+py*ct.vv[j].y+pz); // (made positive for in front of camera)
                if(Tdist<mn*0.995) Tdist = mn; // check here this is buggy
                if(Tdist>mx*1.005) Tdist = mx; // check here this is buggy
                if(Tdist<=0.0) continue;
                xx = static_cast<int>(ct.vv[j].x)-1;
                yy = static_cast<int>(ct.vv[j].y)-1;
                if(xx>=0 && xx<w && yy>=0 && yy<h)
                {
                    if(depthmap[xx*h+yy]>Tdist)
                    {
                        depthmap[xx*h+yy] = Tdist;
                    }
                }
            }
        }
    }

    return 0;
}

/** \brief Orthographic-renders a given mesh to a depthmap and optionally scalarmap
 *
 * \param[in] m Input mesh (must contain faces, face-scalars are preferred)
 * \param[out] depthmap Column-major depth map
 * \param[out] scalarmap Column-major scalar map
 * \param[in] g Camera focal-length (in px)
 * \param[in] h Map height (in px)
 * \param[in] w Map width (in px)
 * \param[in] doscalar Render scalar map (true/false)
 * \param[in] cutoff Render front z-cut-off (default - 0.0)
 * \param[in] eps Render tolerance (default - 0.0001)
 * \return Error code (0-success)
 *
 */

SOFTRENDERER_API int sr_mesh_render_ortho(MESH m, std::vector<FLOATDATA>& depthmap, std::vector<FLOATDATA>& scalarmap, FLOATDATA g, INTDATA h, INTDATA w, bool doscalar, FLOATDATA cutoff, FLOATDATA eps)
{
    INTDATA nv, nf, xx, yy;
    std::vector<FLOATDATA> x, y;
    FLOATDATA Tdist, fw, fh, a, b, c, d, px, py, pd; //, mn, mx;

    fw = w/2.0;
    fh = h/2.0;

    MESH_VERTEX mv = m->vertices;
    nv = m->num_vertices;
    MESH_SCALAR mvs = m->vscalars;
    MESH_FACE mf = m->faces;
    nf = m->num_faces;
    MESH_SCALAR mfs = m->fscalars;

    x.resize(nv);
    y.resize(nv);

    for(INTDATA i=0; i<nv; ++i)
    {
        const MESH_VERTEX mvc = mv+i;
        if(mvc->z<-cutoff)
        {
            x[i] = fw+g*mvc->x;
            y[i] = fh-g*mvc->y;
            //if(x[i]<0.0||y[i]<0.0||x[i]>=w||y[i]>=h) // this is valid only when sufficiently sampled
            if(x[i]<-w/4||y[i]<-h/4||x[i]>=(5*w)/4||y[i]>=(5*h)/4)
            {
                x[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
                y[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
            }
        }
        else
        {
            x[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
            y[i] = std::numeric_limits<FLOATDATA>::quiet_NaN();
        }
    }

    if(m->is_fscalars&&doscalar)
    {
        depthmap.resize(h*w);
        std::fill(depthmap.begin(), depthmap.end(), std::numeric_limits<FLOATDATA>::max());
        scalarmap.resize(h*w);
        std::fill(scalarmap.begin(), scalarmap.end(), std::numeric_limits<FLOATDATA>::quiet_NaN());

        for(INTDATA i=0; i<nf; ++i)
        {
            const INTDATA* mfcv = mf[i].vertices;
            const INTDATA vx1 = mfcv[0], vx2 = mfcv[1], vx3 = mfcv[2];
            if(isnan(x[vx1])||isnan(x[vx2])||isnan(x[vx3])) continue;
            Triangle<FLOATDATA> ct(x[vx1],y[vx1],x[vx2],y[vx2],x[vx3],y[vx3]);
            ct.check();
            ct.draw_triangle();
            a = (mv[vx2].y-mv[vx1].y)*(mv[vx3].z-mv[vx1].z)
                -(mv[vx3].y-mv[vx1].y)*(mv[vx2].z-mv[vx1].z);
            b = (mv[vx2].z-mv[vx1].z)*(mv[vx3].x-mv[vx1].x)
                -(mv[vx3].z-mv[vx1].z)*(mv[vx2].x-mv[vx1].x);
            c = (mv[vx2].x-mv[vx1].x)*(mv[vx3].y-mv[vx1].y)
                -(mv[vx3].x-mv[vx1].x)*(mv[vx2].y-mv[vx1].y);
            d = -(a*mv[vx1].x+b*mv[vx1].y+c*mv[vx1].z);
            px = a/g;
            py = -b/g;
            pd = b*fh/g-a*fw/g+d;
            FLOATDATA mn = -std::max(std::max(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            FLOATDATA mx = -std::min(std::min(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            for(size_t j=0; j<ct.vv.size(); ++j)
            {
                Tdist = (px*ct.vv[j].x+py*ct.vv[j].y+pd)/c; // (made positive for in front of camera)
                if(Tdist<mn*0.995) Tdist = mn; // check here this is buggy
                if(Tdist>mx*1.005) Tdist = mx; // check here this is buggy
                if(Tdist<=0.0) continue;
                xx = static_cast<int>(ct.vv[j].x)-1;
                yy = static_cast<int>(ct.vv[j].y)-1;
                if(xx>=0 && xx<w && yy>=0 && yy<h)
                {
                    if (depthmap[xx*h + yy] > Tdist)
                    {
                        depthmap[xx*h + yy] = Tdist;
                        scalarmap[xx*h + yy] = mfs[i];// assumed faces have scalars
                    }
                }
            }
        }
    }
    if(m->is_vscalars&&doscalar)
    {
        depthmap.resize(h*w);
        std::fill(depthmap.begin(), depthmap.end(), std::numeric_limits<FLOATDATA>::max());
        scalarmap.resize(h*w);
        std::fill(scalarmap.begin(), scalarmap.end(), std::numeric_limits<FLOATDATA>::quiet_NaN());

        for(INTDATA i=0; i<nf; ++i)
        {
            const INTDATA* mfcv = mf[i].vertices;
            const INTDATA vx1 = mfcv[0], vx2 = mfcv[1], vx3 = mfcv[2];
            if(isnan(x[vx1])||isnan(x[vx2])||isnan(x[vx3])) continue;
            Triangle<FLOATDATA> ct(x[vx1],y[vx1],x[vx2],y[vx2],x[vx3],y[vx3]);
            ct.check();
            ct.draw_triangle();
            a = (mv[vx2].y-mv[vx1].y)*(mv[vx3].z-mv[vx1].z)
                -(mv[vx3].y-mv[vx1].y)*(mv[vx2].z-mv[vx1].z);
            b = (mv[vx2].z-mv[vx1].z)*(mv[vx3].x-mv[vx1].x)
                -(mv[vx3].z-mv[vx1].z)*(mv[vx2].x-mv[vx1].x);
            c = (mv[vx2].x-mv[vx1].x)*(mv[vx3].y-mv[vx1].y)
                -(mv[vx3].x-mv[vx1].x)*(mv[vx2].y-mv[vx1].y);
            d = -(a*mv[vx1].x+b*mv[vx1].y+c*mv[vx1].z);
            px = a/g;
            py = -b/g;
            pd = b*fh/g-a*fw/g+d;
            FLOATDATA mn = -std::max(std::max(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            FLOATDATA mx = -std::min(std::min(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            for(size_t j=0; j<ct.vv.size(); ++j)
            {
                Tdist = (px*ct.vv[j].x+py*ct.vv[j].y+pd)/c; // (made positive for in front of camera)
                if(Tdist<mn*0.995) Tdist = mn; // check here this is buggy
                if(Tdist>mx*1.005) Tdist = mx; // check here this is buggy
                if(Tdist<=0.0) continue;
                xx = static_cast<int>(ct.vv[j].x)-1;
                yy = static_cast<int>(ct.vv[j].y)-1;
                if(xx>=0 && xx<w && yy>=0 && yy<h)
                {
                    if (depthmap[xx*h + yy] > Tdist)
                    {
                        depthmap[xx*h + yy] = Tdist;
                        scalarmap[xx*h + yy] = mvs[vx1];// assumed all vertex has same scalar
                    }
                }
            }
        }
    }
    else
    {
        depthmap.resize(h*w);
        std::fill(depthmap.begin(), depthmap.end(), std::numeric_limits<FLOATDATA>::max());

        for(INTDATA i=0; i<nf; ++i)
        {
            const INTDATA* mfcv = mf[i].vertices;
            const INTDATA vx1 = mfcv[0], vx2 = mfcv[1], vx3 = mfcv[2];
            if(isnan(x[vx1])||isnan(x[vx2])||isnan(x[vx3])) continue;
            Triangle<FLOATDATA> ct(x[vx1],y[vx1],x[vx2],y[vx2],x[vx3],y[vx3]);
            ct.check();
            ct.draw_triangle();
            a = (mv[vx2].y-mv[vx1].y)*(mv[vx3].z-mv[vx1].z)
                -(mv[vx3].y-mv[vx1].y)*(mv[vx2].z-mv[vx1].z);
            b = (mv[vx2].z-mv[vx1].z)*(mv[vx3].x-mv[vx1].x)
                -(mv[vx3].z-mv[vx1].z)*(mv[vx2].x-mv[vx1].x);
            c = (mv[vx2].x-mv[vx1].x)*(mv[vx3].y-mv[vx1].y)
                -(mv[vx3].x-mv[vx1].x)*(mv[vx2].y-mv[vx1].y);
            d = -(a*mv[vx1].x+b*mv[vx1].y+c*mv[vx1].z);
            px = a/g;
            py = -b/g;
            pd = b*fh/g-a*fw/g+d;
            FLOATDATA mn = -std::max(std::max(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            FLOATDATA mx = -std::min(std::min(mv[vx1].z,mv[vx2].z),mv[vx3].z);
            for(size_t j=0; j<ct.vv.size(); ++j)
            {
                Tdist = (px*ct.vv[j].x+py*ct.vv[j].y+pd)/c; // (made positive for in front of camera)
                if(Tdist<mn*0.995) Tdist = mn; // check here this is buggy
                if(Tdist>mx*1.005) Tdist = mx; // check here this is buggy
                if(Tdist<=0.0) continue;
                xx = static_cast<int>(ct.vv[j].x)-1;
                yy = static_cast<int>(ct.vv[j].y)-1;
                if(xx>=0 && xx<w && yy>=0 && yy<h)
                {
                    if(depthmap[xx*h+yy]>Tdist)
                    {
                        depthmap[xx*h+yy] = Tdist;
                    }
                }
            }
        }
    }

    return 0;
}

