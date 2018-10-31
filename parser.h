//FATIMA REY
#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <iostream>

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    

    struct Vec3f
    {
        float x, y, z;
        Vec3f() : x(0.0f), y(0.0f), z(0.0f) {}
        Vec3f(float X, float Y, float Z) : x(X), y(Y), z(Z) {}

        //basic vector operations
        Vec3f operator + (Vec3f v2) {return Vec3f(x+v2.x, y+v2.y, z+v2.z); }
        Vec3f operator - (Vec3f v2) {return Vec3f(x-v2.x, y-v2.y, z-v2.z); }
        Vec3f operator * (float k) {return Vec3f(x*k, y*k, z*k); }
        Vec3f operator / (float k) {return Vec3f(x/k, y/k, z/k); }

        void operator = (Vec3f v2) {x=v2.x; y=v2.y; z=v2.z;}
        
        bool operator == (Vec3f v2) { return x==v2.x && y==v2.y && z== v2.z ;} 

        float dot(Vec3f v2) {return (x*v2.x + y*v2.y + z*v2.z);}

        Vec3f cross(Vec3f w)
        {
            Vec3f u;
            u.x = y*w.z - w.y*z;
            u.y = w.x*z - x*w.z;
            u.z = x*w.y - w.x*y;
            return u;
        }
        
        float magnitude() { return sqrtf(x*x + y*y +z *z); }

        Vec3f normalize()
        {
            float magnitude = sqrtf(x*x + y*y + z*z);
            return Vec3f(x/magnitude, y/magnitude, z/magnitude);
        }

    };

    struct Ray 
    {
        Vec3f origin, direction;
        Ray() : origin(0.0f,0.0f,0.0f), direction(0.0f,0.0f,0.0f) {}
        Ray(Vec3f o, Vec3f d) {origin=o; direction=d;}

        void operator = (Ray r2) {origin=r2.origin; direction=r2.direction; }
    };

    struct Matrix
    {
        float a11, a12, a13, a14;       //first row
        float a21, a22, a23, a24;       //second row
        float a31, a32, a33, a34;       //third row
        float a41, a42, a43, a44;       //fourth row

        Matrix()
        {
            a11=0; a12=0; a13=0; a14=0;     //first row
            a21=0; a22=0; a23=0; a24=0;     //second row
            a31=0; a32=0; a33=0; a34=0;     //third row
            a41=0; a42=0; a43=0; a44=0;
        }

        Matrix( float c11, float c12, float c13, float c14,
                float c21, float c22, float c23, float c24,
                float c31, float c32, float c33, float c34,
                float c41, float c42, float c43, float c44) 
        {
            a11=c11; a12=c12; a13=c13; a14=c14;     //first row
            a21=c21; a22=c22; a23=c23; a24=c24;     //second row
            a31=c31; a32=c32; a33=c33; a34=c34;     //third row
            a41=c41; a42=c42; a43=c43; a44=c44;
        }
    
    
        Matrix operator * (Matrix m2) 
        {
            float c11, c12, c13, c14;
            float c21, c22, c23, c24;
            float c31, c32, c33, c34;
            float c41, c42, c43, c44;

            c11 = a11*m2.a11 + a12*m2.a21 + a13*m2.a31 + a14*m2.a41;
            c12 = a11*m2.a12 + a12*m2.a22 + a13*m2.a32 + a14*m2.a42;
            c13 = a11*m2.a13 + a12*m2.a23 + a13*m2.a33 + a14*m2.a43;
            c14 = a11*m2.a14 + a12*m2.a24 + a13*m2.a34 + a14*m2.a44;

            c21 = a21*m2.a11 + a22*m2.a21 + a23*m2.a31 + a24*m2.a41;
            c22 = a21*m2.a12 + a22*m2.a22 + a23*m2.a32 + a24*m2.a42;
            c23 = a21*m2.a13 + a22*m2.a23 + a23*m2.a33 + a24*m2.a43;
            c24 = a21*m2.a14 + a22*m2.a24 + a23*m2.a34 + a24*m2.a44;

            c31 = a31*m2.a11 + a32*m2.a21 + a33*m2.a31 + a34*m2.a41;
            c32 = a31*m2.a12 + a32*m2.a22 + a33*m2.a32 + a34*m2.a42;
            c33 = a31*m2.a13 + a32*m2.a23 + a33*m2.a33 + a34*m2.a43;
            c34 = a31*m2.a14 + a32*m2.a24 + a33*m2.a34 + a34*m2.a44;

            c41 = a41*m2.a11 + a42*m2.a21 + a43*m2.a31 + a44*m2.a41;
            c42 = a41*m2.a12 + a42*m2.a22 + a43*m2.a32 + a44*m2.a42;
            c43 = a41*m2.a13 + a42*m2.a23 + a43*m2.a33 + a44*m2.a43;
            c44 = a41*m2.a14 + a42*m2.a24 + a43*m2.a34 + a44*m2.a44;
        
            return Matrix(c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,c41,c42,c43,c44); }

        Vec3f operator * (Vec3f v) 
        {
            float rx, ry, rz;

            rx = a11*v.x + a12*v.y + a13*v.z + a14;
            ry = a21*v.x + a22*v.y + a23*v.z + a24;
            rz = a31*v.x + a32*v.y + a33*v.z + a34;
        
            return Vec3f(rx, ry, rz);
        }

        Vec3f vecMultip(Vec3f v)
        {
            float rx, ry, rz;

            rx = a11*v.x + a12*v.y + a13*v.z ;
            ry = a21*v.x + a22*v.y + a23*v.z ;
            rz = a31*v.x + a32*v.y + a33*v.z ;
        
            return Vec3f(rx, ry, rz);
        }
    
        float determinant()
        {
            float d = a11*a22*a33*a44 + a11*a23*a34*a42 + a11*a24*a32*a43 
                    + a12*a21*a34*a43 + a12*a23*a31*a44 + a12*a24*a33*a41 
                    + a13*a21*a32*a44 + a13*a22*a34*a41 + a13*a24*a31*a42
                    + a14*a21*a33*a42 + a14*a22*a31*a43 + a14*a23*a32*a41
                    - a11*a22*a34*a43 - a11*a23*a32*a44 - a11*a24*a33*a42
                    - a12*a21*a33*a44 - a12*a23*a34*a41 - a12*a24*a31*a43
                    - a13*a21*a34*a42 - a13*a22*a31*a44 - a13*a24*a32*a41
                    - a14*a21*a32*a43 - a14*a22*a33*a41 - a14*a23*a31*a42;

            return d;
        }

        Matrix inverse()
        {
            float c11, c12, c13, c14;
            float c21, c22, c23, c24;
            float c31, c32, c33, c34;
            float c41, c42, c43, c44;

            c11 = a22*a33*a44 + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*a44 - a24*a33*a42; 
            c12 = a12*a34*a43 + a13*a32*a44 + a14*a33*a42 - a12*a33*a44 - a13*a34*a42 - a14*a32*a43;
            c13 = a12*a23*a44 + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*a44 - a14*a23*a42;
            c14 = a12*a24*a33 + a13*a22*a34 + a14*a23*a32 - a12*a23*a34 - a13*a24*a32 - a14*a22*a33;

            c21 = a21*a34*a43 + a23*a31*a44 + a24*a33*a41 - a21*a33*a44 - a23*a34*a41 - a24*a31*a43;
            c22 = a11*a33*a44 + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*a44 - a14*a33*a41;
            c23 = a11*a24*a43 + a13*a21*a44 + a14*a23*a41 - a11*a23*a44 - a13*a24*a41 - a14*a21*a43;
            c24 = a11*a23*a34 + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*a34 - a14*a23*a31;

            c31 = a21*a32*a44 + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*a44 - a24*a32*a41;
            c32 = a11*a34*a42 + a12*a31*a44 + a14*a32*a41 - a11*a32*a44 - a12*a34*a41 - a14*a31*a42;
            c33 = a11*a22*a44 + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*a44 - a14*a22*a41;
            c34 = a11*a24*a32 + a12*a21*a34 + a14*a22*a31 - a11*a22*a34 - a12*a24*a31 - a14*a21*a32;

            c41 = a21*a33*a42 + a22*a31*a43 + a23*a32*a41 - a21*a32*a43 - a22*a33*a41 - a23*a31*a42;
            c42 = a11*a32*a43 + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*a43 - a13*a32*a41;
            c43 = a11*a23*a42 + a12*a21*a43 + a13*a22*a41 - a11*a22*a43 - a12*a23*a41 - a13*a21*a42;
            c44 = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*a33 - a13*a22*a31;

            float d = this->determinant();
        
            d = 1/d;
        
            c11 *= d; c12 *= d; c13 *= d; c14 *= d;
            c21 *= d; c22 *= d; c23 *= d; c24 *= d;
            c31 *= d; c32 *= d; c33 *= d; c34 *= d;
            c41 *= d; c42 *= d; c43 *= d; c44 *= d;

            return Matrix(c11,c12,c13,c14,c21,c22,c23,c24,c31,c32,c33,c34,c41,c42,c43,c44);
        }

        void operator = (Matrix m2)
        {

            a11 = m2.a11; a12 = m2.a12; a13=m2.a13; a14=m2.a14;
            a21 = m2.a21; a22 = m2.a22; a23=m2.a23; a24=m2.a24;
            a31 = m2.a31; a32 = m2.a32; a33=m2.a33; a34=m2.a34;
            a41 = m2.a41; a42 = m2.a42; a43=m2.a43; a44=m2.a44;
        }

        Matrix transpose()
        {
            return Matrix(a11,a21,a31,a41,
                          a12,a22,a32,a42,
                          a13,a23,a33,a43,
                          a14,a24,a34,a44);
        }


        Ray operator * (Ray r)
        {
            Vec3f o, d;

            o.x = a11*r.origin.x + a12*r.origin.y + a13*r.origin.z + a14;
            o.y = a21*r.origin.x + a22*r.origin.y + a23*r.origin.z + a24;
            o.z = a31*r.origin.x + a32*r.origin.y + a33*r.origin.z + a34;

            d.x = a11*r.direction.x + a12*r.direction.y + a13*r.direction.z;
            d.y = a21*r.direction.x + a22*r.direction.y + a23*r.direction.z;
            d.z = a31*r.direction.x + a32*r.direction.y + a33*r.direction.z;
           
            return Ray(o, d);
        }

        void print()
        {
            std::cout << a11 << " " << a12 << " " << a13 << " " << a14 << std::endl;
            std::cout << a21 << " " << a22 << " " << a23 << " " << a24 << std::endl;
            std::cout << a31 << " " << a32 << " " << a33 << " " << a34 << std::endl;
            std::cout << a41 << " " << a42 << " " << a43 << " " << a44 << std::endl;
            std::cout << std::endl;
        }

    };


    

    struct Vec3i
    {
        int x, y, z;
    };

    struct Vec4f
    {
        float x, y, z, w;
    };

    struct Camera
    {
        Vec3f position;
        Vec3f gaze;
        Vec3f up;
        Vec4f near_plane;
        float near_distance;
        int image_width, image_height;
        std::string image_name;
    };

    struct PointLight
    {
        Vec3f position;
        Vec3f intensity;
    };

    struct Material
    {
        Vec3f ambient;
        Vec3f diffuse;
        Vec3f specular;
        Vec3f mirror;
        float phong_exponent;
    };

    struct Face
    {
        int v0_id;
        int v1_id;
        int v2_id;

        //khsfkjhaksd
        Face(): v0_id(0),  v1_id(0), v2_id(0) {}
        Face(int a, int b, int c): v0_id(a),  v1_id(b), v2_id(c) {}

        void operator = (Face f) { v0_id=f.v0_id;  v1_id=f.v1_id; v2_id=f.v2_id; }
    };

struct Triangle
    {
        int material_id;
        Face indices;

        Vec3f vertex0, vertex1, vertex2 ;
        Vec3f normal;
        Vec3f unit_normal;

        Material mat;

        Triangle(): material_id(0) , indices(Face(0,0,0)) {}
        Triangle(int m, Face i): material_id(m) , indices(i) {}

        void compute_normal()
        {
            normal = (vertex1 - vertex0).cross(vertex2 - vertex0);
            unit_normal = normal.normalize() ;
        }
 
        float compute_det(float a,float b,float c,float d,float e,float f,float g,float h,float i)
        {
            float det = a*(e*i-h*f) + b*(g*f-d*i) + c*(d*h-e*g);
            return  det ;
        }

        bool is_intersect(Ray& ray, float& t)
        {
            // check if point is in triangle's plane
            if(unit_normal.dot(ray.direction) == 0) 
            {
                return false;
            }

            float detA     = compute_det(vertex0.x-vertex1.x, vertex0.x-vertex2.x, ray.direction.x,
                                         vertex0.y-vertex1.y, vertex0.y-vertex2.y, ray.direction.y,
                                         vertex0.z-vertex1.z, vertex0.z-vertex2.z, ray.direction.z);

            float detalpha = compute_det(vertex0.x-ray.origin.x, vertex0.x-vertex2.x, ray.direction.x,
                                         vertex0.y-ray.origin.y, vertex0.y-vertex2.y, ray.direction.y, 
                                         vertex0.z-ray.origin.z, vertex0.z-vertex2.z, ray.direction.z);

            float detbeta = compute_det(vertex0.x-vertex1.x, vertex0.x-ray.origin.x, ray.direction.x, 
                                        vertex0.y-vertex1.y, vertex0.y-ray.origin.y, ray.direction.y,
                                        vertex0.z-vertex1.z, vertex0.z-ray.origin.z, ray.direction.z);

            float detT     = compute_det(vertex0.x-vertex1.x, vertex0.x-vertex2.x, vertex0.x-ray.origin.x,
                                         vertex0.y-vertex1.y, vertex0.y-vertex2.y, vertex0.y-ray.origin.y,
                                         vertex0.z-vertex1.z, vertex0.z-vertex2.z, vertex0.z-ray.origin.z);

            float alpha = detalpha / detA;
            float beta  = detbeta / detA;

            t = detT / detA;
            if (t<=0) return 0;

            if(alpha >= 0.0f && beta >= 0.0f && alpha+beta <= 1.0f)
                return true;
            else
                return false;
        }
    };
    
    struct Mesh
    {
        int material_id;
        std::vector<Face> faces;

        std::vector<Triangle> mtriangles;




    };

    

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;


        Vec3f center;

        Vec3f unit_normal;

        Material mat;
        
        void compute_normal(Vec3f a)
        {
            unit_normal = (a-center)/radius;
        }   

    

        bool is_intersect(Ray& ray, float& t)
        {
            Vec3f o = ray.origin;
            Vec3f d = ray.direction;

                  

            Vec3f c = this->center;
            float r = this->radius;

            Vec3f oc = o-c;

            float B = 2*(d.dot(oc));
            float C = oc.dot(oc) - r*r ;

            float delta = B*B - 4*C;

            if(delta < -0.00001f)
            {   return false;}

            else
            {               
                float t1 = (-B - sqrt(delta)) / 2;
                float t2 = (-B + sqrt(delta)) / 2;

                t = (t1 < t2) ? t1 : t2 ;

                if (t<=0) return 0; 

                Vec3f p = o + d*t ;

                compute_normal(p);

               

                return true;    
            }

        }


    };

    struct Scene
    {
        //Data
        Vec3i background_color;
        float shadow_ray_epsilon;
        int max_recursion_depth;
        std::vector<Camera> cameras;
        Vec3f ambient_light;
        std::vector<PointLight> point_lights;
        std::vector<Material> materials;
        std::vector<Vec3f> vertex_data;
        std::vector<Mesh> meshes;
        std::vector<Triangle> triangles;
        std::vector<Sphere> spheres;

        //Functions
        void loadFromXml(const std::string& filepath);
		bool isIntersected(Ray ray, float& t, Material& imat, Vec3f& un);
		Vec3i computeAmbientLight(Ray ray, float& t, Material& material, Vec3f& un, int& count);

    };

}

#endif
