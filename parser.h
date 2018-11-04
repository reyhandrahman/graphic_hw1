//FATIMA REY
#ifndef __HW1__PARSER__
#define __HW1__PARSER__

#include <string>
#include <vector>
#include <limits>
<<<<<<< HEAD
=======
#include <cmath>
#include <iostream>
>>>>>>> 7bafb10aa13de0f7c8b9cf833df1e33240e7eace

namespace parser
{
    //Notice that all the structures are as simple as possible
    //so that you are not enforced to adopt any style or design.
    

    struct Vec3f
    {
        float x, y, z;

        //constructors 
        Vec3f() : x(0.0f), y(0.0f), z(0.0f) {}
        Vec3f(float X, float Y, float Z) : x(X), y(Y), z(Z) {}

        //new operators to help
        Vec3f operator + (Vec3f v2) {return Vec3f(x+v2.x, y+v2.y, z+v2.z); }
        Vec3f operator - (Vec3f v2) {return Vec3f(x-v2.x, y-v2.y, z-v2.z); }
        Vec3f operator * (float k) {return Vec3f(x*k, y*k, z*k); }
        Vec3f operator / (float k) {return Vec3f(x/k, y/k, z/k); }

        void operator = (Vec3f v2) {x=v2.x; y=v2.y; z=v2.z;}
        bool operator == (Vec3f v2) { return x==v2.x && y==v2.y && z== v2.z ;} 

        //dot and cross products
        float dot(Vec3f v2) 
        {
            return (x*v2.x + y*v2.y + z*v2.z);
        }

        Vec3f cross(Vec3f w)
        {
            Vec3f u;
            u.x = y*w.z - w.y*z;
            u.y = w.x*z - x*w.z;
            u.z = x*w.y - w.x*y;
            return u;
        }

        Vec3f normalize()
        {
            float magnitude = sqrtf(x*x + y*y + z*z);
            return Vec3f(x/magnitude, y/magnitude, z/magnitude);
        }
    };

    struct Ray 
    {
        Vec3f o, d;
        //empty ray
        Ray() : o(0.0f,0.0f,0.0f), d(0.0f,0.0f,0.0f) {}

        //defined ray
        Ray(Vec3f orgn, Vec3f drct) {o=orgn; d=drct;}

        void operator = (Ray r2) {o=r2.o; d=r2.d; }
    };

    

    struct Matrix
    {
        float intersectionArray[4][4];
        
        Matrix()
        {
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    intersectionArray[i][j]=0;
                }
            }
        }
           
            Matrix(float a11,float a12,float a13, float a14,       
            float a21,float a22,float a23, float a24,       
            float a31, float a32, float a33,float a34,       
            float a41, float a42,float a43, float a44)
            {
               intersectionArray[0][0]=a11;
               intersectionArray[0][1]=a12;
               intersectionArray[0][2]=a13;
               intersectionArray[0][3]=a14;
               intersectionArray[1][0]=a21;
               intersectionArray[1][1]=a22;
               intersectionArray[1][2]=a23;
               intersectionArray[1][3]=a24;
               intersectionArray[2][0]=a31;
               intersectionArray[2][1]=a32;
               intersectionArray[2][2]=a33;
               intersectionArray[2][3]=a34;
               intersectionArray[3][0]=a41;
               intersectionArray[3][1]=a42;
               intersectionArray[3][2]=a43;
               intersectionArray[3][3]=a44;
            }

          Ray operator * (Ray r)
        {
            Vec3f o, d;

            o.x = intersectionArray[0][0]*r.o.x + intersectionArray[0][1]*r.o.y + intersectionArray[0][2]*r.o.z + intersectionArray[0][3];
            o.y = intersectionArray[1][0]*r.o.x + intersectionArray[1][1]*r.o.y + intersectionArray[1][2]*r.o.z + intersectionArray[1][3];
            o.z = intersectionArray[2][0]*r.o.x + intersectionArray[2][1]*r.o.y + intersectionArray[2][2]*r.o.z + intersectionArray[2][3];

            d.x = intersectionArray[0][0]*r.d.x + intersectionArray[0][1]*r.d.y + intersectionArray[0][2]*r.d.z;
            d.y = intersectionArray[1][0]*r.d.x + intersectionArray[1][1]*r.d.y + intersectionArray[1][2]*r.d.z;
            d.z = intersectionArray[2][0]*r.d.x + intersectionArray[2][1]*r.d.y + intersectionArray[2][2]*r.d.z;

            return Ray(o, d);
        }


         Matrix operator * (Matrix m2) 
        {
            
            Matrix result;
          
            result.intersectionArray[0][0] = intersectionArray[0][0]*m2.intersectionArray[0][0] + intersectionArray[0][1]*m2.intersectionArray[1][0] + intersectionArray[0][2]*m2.intersectionArray[2][0] + intersectionArray[0][3]*m2.intersectionArray[3][0];
            result.intersectionArray[0][1] = intersectionArray[0][0]*m2.intersectionArray[0][1] + intersectionArray[0][1]*m2.intersectionArray[1][1] + intersectionArray[0][2]*m2.intersectionArray[2][1] + intersectionArray[0][3]*m2.intersectionArray[3][1];
            result.intersectionArray[0][2]= intersectionArray[0][0]*m2.intersectionArray[0][2] + intersectionArray[0][1]*m2.intersectionArray[1][2] + intersectionArray[0][2]*m2.intersectionArray[2][2] + intersectionArray[0][3]*m2.intersectionArray[3][2];
            result.intersectionArray[0][3]= intersectionArray[0][0]*m2.intersectionArray[0][3] + intersectionArray[0][1]*m2.intersectionArray[1][3] + intersectionArray[0][2]*m2.intersectionArray[2][3] + intersectionArray[0][3]*m2.intersectionArray[3][3];

            result.intersectionArray[1][0] = intersectionArray[1][0]*m2.intersectionArray[0][0] + intersectionArray[1][1]*m2.intersectionArray[1][0] + intersectionArray[1][2]*m2.intersectionArray[2][0] + intersectionArray[1][3]*m2.intersectionArray[3][0];
            result.intersectionArray[1][1] = intersectionArray[1][0]*m2.intersectionArray[0][1] + intersectionArray[1][1]*m2.intersectionArray[1][1] + intersectionArray[1][2]*m2.intersectionArray[2][1] + intersectionArray[1][3]*m2.intersectionArray[3][1];
            result.intersectionArray[1][2]= intersectionArray[1][0]*m2.intersectionArray[0][2] + intersectionArray[1][1]*m2.intersectionArray[1][2] + intersectionArray[1][2]*m2.intersectionArray[2][2] + intersectionArray[1][3]*m2.intersectionArray[3][2];
            result.intersectionArray[1][3]= intersectionArray[1][0]*m2.intersectionArray[0][3] + intersectionArray[1][1]*m2.intersectionArray[1][3] + intersectionArray[1][2]*m2.intersectionArray[2][3] + intersectionArray[1][3]*m2.intersectionArray[3][3];

            result.intersectionArray[2][0] = intersectionArray[2][0]*m2.intersectionArray[0][0] + intersectionArray[2][1]*m2.intersectionArray[1][0] + intersectionArray[2][2]*m2.intersectionArray[2][0] + intersectionArray[2][3]*m2.intersectionArray[3][0];
            result.intersectionArray[2][1] = intersectionArray[2][0]*m2.intersectionArray[0][1] + intersectionArray[2][1]*m2.intersectionArray[1][1] + intersectionArray[2][2]*m2.intersectionArray[2][1] + intersectionArray[2][3]*m2.intersectionArray[3][1];
            result.intersectionArray[2][2]= intersectionArray[2][0]*m2.intersectionArray[0][2] + intersectionArray[2][1]*m2.intersectionArray[1][2] + intersectionArray[2][2]*m2.intersectionArray[2][2] + intersectionArray[2][3]*m2.intersectionArray[3][2];
            result.intersectionArray[2][3]= intersectionArray[2][0]*m2.intersectionArray[0][3] + intersectionArray[2][1]*m2.intersectionArray[1][3] + intersectionArray[2][2]*m2.intersectionArray[2][3] + intersectionArray[2][3]*m2.intersectionArray[3][3];

            result.intersectionArray[3][0] = intersectionArray[3][0]*m2.intersectionArray[0][0] + intersectionArray[3][1]*m2.intersectionArray[1][0] + intersectionArray[3][2]*m2.intersectionArray[2][0] + intersectionArray[3][3]*m2.intersectionArray[3][0];
            result.intersectionArray[3][1] = intersectionArray[3][0]*m2.intersectionArray[0][1] + intersectionArray[3][1]*m2.intersectionArray[1][1] + intersectionArray[3][2]*m2.intersectionArray[2][1] + intersectionArray[3][3]*m2.intersectionArray[3][1];
            result.intersectionArray[3][2]= intersectionArray[3][0]*m2.intersectionArray[0][2] + intersectionArray[3][1]*m2.intersectionArray[1][2] + intersectionArray[3][2]*m2.intersectionArray[2][2] + intersectionArray[3][3]*m2.intersectionArray[3][2];
            result.intersectionArray[3][3]= intersectionArray[3][0]*m2.intersectionArray[0][3] + intersectionArray[3][1]*m2.intersectionArray[1][3] + intersectionArray[3][2]*m2.intersectionArray[2][3] + intersectionArray[3][3]*m2.intersectionArray[3][3];

            return result;
        }
        

        Vec3f operator * (Vec3f v) 
        {
            float rx, ry, rz;

            rx = this->intersectionArray[0][0]*v.x + this->intersectionArray[0][1]*v.y + this->intersectionArray[0][2]*v.z + this->intersectionArray[0][3];
            ry = this->intersectionArray[1][0]*v.x + this->intersectionArray[1][1]*v.y + this->intersectionArray[1][2]*v.z + this->intersectionArray[1][3];
            rz = this->intersectionArray[2][0]*v.x + this->intersectionArray[2][1]*v.y + this->intersectionArray[2][2]*v.z + this->intersectionArray[2][3];
        
            return Vec3f(rx, ry, rz);
        }

        Vec3f vecMultip(Vec3f v)
        {
            float rx, ry, rz;

           rx = this->intersectionArray[0][0]*v.x + this->intersectionArray[0][1]*v.y + this->intersectionArray[0][2]*v.z;
            ry = this->intersectionArray[1][0]*v.x + this->intersectionArray[1][1]*v.y + this->intersectionArray[1][2]*v.z;
            rz = this->intersectionArray[2][0]*v.x + this->intersectionArray[2][1]*v.y + this->intersectionArray[2][2]*v.z;
        
            return Vec3f(rx, ry, rz);
        }
            

        float determinant()
        {
            float det = intersectionArray[0][0]*intersectionArray[1][1]*intersectionArray[2][2]*intersectionArray[3][3]+ intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[2][3]*intersectionArray[3][1]+intersectionArray[0][0]*intersectionArray[1][3]*intersectionArray[2][1]*intersectionArray[3][2]
                    +intersectionArray[0][1]*intersectionArray[1][0]*intersectionArray[2][3]*intersectionArray[3][2]+ intersectionArray[0][1]*intersectionArray[1][2]*intersectionArray[2][0]*intersectionArray[3][3]+intersectionArray[0][1]*intersectionArray[1][3]*intersectionArray[2][2]*intersectionArray[3][0]
                    +intersectionArray[0][2]*intersectionArray[1][0]*intersectionArray[2][1]*intersectionArray[3][3]+ intersectionArray[0][2]*intersectionArray[1][1]*intersectionArray[2][3]*intersectionArray[3][0]+intersectionArray[0][2]*intersectionArray[1][3]*intersectionArray[2][0]*intersectionArray[3][1]
                    +intersectionArray[0][3]*intersectionArray[1][0]*intersectionArray[2][2]*intersectionArray[3][1]+ intersectionArray[0][3]*intersectionArray[1][1]*intersectionArray[2][0]*intersectionArray[3][2]+intersectionArray[0][3]*intersectionArray[1][2]*intersectionArray[2][1]*intersectionArray[3][0]
                    -intersectionArray[0][0]*intersectionArray[1][1]*intersectionArray[2][3]*intersectionArray[3][2]- intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[2][1]*intersectionArray[3][3]-intersectionArray[0][0]*intersectionArray[1][3]*intersectionArray[2][2]*intersectionArray[3][1]                    
                    -intersectionArray[0][1]*intersectionArray[1][0]*intersectionArray[2][2]*intersectionArray[3][3]- intersectionArray[0][1]*intersectionArray[1][2]*intersectionArray[2][3]*intersectionArray[3][0]-intersectionArray[0][1]*intersectionArray[1][3]*intersectionArray[2][0]*intersectionArray[3][2]
                    -intersectionArray[0][2]*intersectionArray[1][0]*intersectionArray[2][3]*intersectionArray[3][1]- intersectionArray[0][2]*intersectionArray[1][1]*intersectionArray[2][0]*intersectionArray[3][3]-intersectionArray[0][2]*intersectionArray[1][3]*intersectionArray[2][1]*intersectionArray[3][0]
                    -intersectionArray[0][3]*intersectionArray[1][0]*intersectionArray[2][1]*intersectionArray[3][2]- intersectionArray[0][3]*intersectionArray[1][1]*intersectionArray[2][2]*intersectionArray[3][0]-intersectionArray[0][3]*intersectionArray[1][2]*intersectionArray[2][0]*intersectionArray[3][1];
            return det;
        }

        Matrix inverse()
        {
            Matrix result;
            result.intersectionArray[0][0] = intersectionArray[1][1]*intersectionArray[2][2]*intersectionArray[3][3]+intersectionArray[1][2]*intersectionArray[2][3]*intersectionArray[3][1]+intersectionArray[1][3]*intersectionArray[2][1]*intersectionArray[3][2]-intersectionArray[1][1]*intersectionArray[2][3]*intersectionArray[3][2]-intersectionArray[1][2]*intersectionArray[2][1]*intersectionArray[3][3]-intersectionArray[1][3]*intersectionArray[2][2]*intersectionArray[3][1];
            result.intersectionArray[0][1] = intersectionArray[0][1]*intersectionArray[2][3]*intersectionArray[3][2]+intersectionArray[0][2]*intersectionArray[2][1]*intersectionArray[3][3]+intersectionArray[0][3]*intersectionArray[2][2]*intersectionArray[3][3]-intersectionArray[0][1]*intersectionArray[2][2]*intersectionArray[3][3]-intersectionArray[0][2]*intersectionArray[2][3]*intersectionArray[3][1]-intersectionArray[0][3]*intersectionArray[2][1]*intersectionArray[3][2];
            result.intersectionArray[0][2]= intersectionArray[0][1]*intersectionArray[1][2]*intersectionArray[3][3]+intersectionArray[0][2]*intersectionArray[1][3]*intersectionArray[3][1]+intersectionArray[0][3]*intersectionArray[1][1]*intersectionArray[3][2]-intersectionArray[0][1]*intersectionArray[1][3]*intersectionArray[3][2]-intersectionArray[0][2]*intersectionArray[1][1]*intersectionArray[3][3]-intersectionArray[0][3]*intersectionArray[1][2]*intersectionArray[3][1];
            result.intersectionArray[0][3]= intersectionArray[0][1]*intersectionArray[1][3]*intersectionArray[2][2]+intersectionArray[0][2]*intersectionArray[1][1]*intersectionArray[2][3]+intersectionArray[0][3]*intersectionArray[1][2]*intersectionArray[2][1]-intersectionArray[0][1]*intersectionArray[1][2]*intersectionArray[2][3]-intersectionArray[0][2]*intersectionArray[1][3]*intersectionArray[2][1]-intersectionArray[0][3]*intersectionArray[1][1]*intersectionArray[2][2];

            result.intersectionArray[1][0] = intersectionArray[1][0]*intersectionArray[2][3]*intersectionArray[3][2]+intersectionArray[1][2]*intersectionArray[2][0]*intersectionArray[3][3]+intersectionArray[1][3]*intersectionArray[2][2]*intersectionArray[3][0]-intersectionArray[1][0]*intersectionArray[2][2]*intersectionArray[3][3]-intersectionArray[1][2]*intersectionArray[2][3]*intersectionArray[3][0]-intersectionArray[1][3]*intersectionArray[2][0]*intersectionArray[3][2];
            result.intersectionArray[1][1] = intersectionArray[0][0]*intersectionArray[2][2]*intersectionArray[3][3]+intersectionArray[0][2]*intersectionArray[2][3]*intersectionArray[3][0]+intersectionArray[0][3]*intersectionArray[2][0]*intersectionArray[3][2]-intersectionArray[0][0]*intersectionArray[2][3]*intersectionArray[3][2]-intersectionArray[0][2]*intersectionArray[2][0]*intersectionArray[3][3]-intersectionArray[0][3]*intersectionArray[2][2]*intersectionArray[3][0];
            result.intersectionArray[1][2]= intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[3][2]+intersectionArray[0][2]*intersectionArray[1][0]*intersectionArray[3][3]+intersectionArray[0][3]*intersectionArray[1][2]*intersectionArray[3][0]-intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[3][3]-intersectionArray[0][2]*intersectionArray[1][3]*intersectionArray[3][0]-intersectionArray[0][3]*intersectionArray[1][0]*intersectionArray[3][2];
            result.intersectionArray[1][3]= intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[2][3]+intersectionArray[0][2]*intersectionArray[1][3]*intersectionArray[2][0]+intersectionArray[0][3]*intersectionArray[1][0]*intersectionArray[2][2]-intersectionArray[0][0]*intersectionArray[1][3]*intersectionArray[2][2]-intersectionArray[0][2]*intersectionArray[1][0]*intersectionArray[3][3]-intersectionArray[0][3]*intersectionArray[1][2]*intersectionArray[2][0];

            result.intersectionArray[2][0] = intersectionArray[1][0]*intersectionArray[2][1]*intersectionArray[3][3]+intersectionArray[1][1]*intersectionArray[2][3]*intersectionArray[3][0]+intersectionArray[1][3]*intersectionArray[2][0]*intersectionArray[3][1]-intersectionArray[1][0]*intersectionArray[2][3]*intersectionArray[3][1]-intersectionArray[1][1]*intersectionArray[2][0]*intersectionArray[3][3]-intersectionArray[1][3]*intersectionArray[2][1]*intersectionArray[3][0];
            result.intersectionArray[2][1] = intersectionArray[0][0]*intersectionArray[2][3]*intersectionArray[3][1]+intersectionArray[0][1]*intersectionArray[2][0]*intersectionArray[3][3]+intersectionArray[0][3]*intersectionArray[2][1]*intersectionArray[3][3]-intersectionArray[0][0]*intersectionArray[2][1]*intersectionArray[3][3]-intersectionArray[0][1]*intersectionArray[2][3]*intersectionArray[3][0]-intersectionArray[0][3]*intersectionArray[2][0]*intersectionArray[3][1];
            result.intersectionArray[2][2]= intersectionArray[0][0]*intersectionArray[1][1]*intersectionArray[3][3]+intersectionArray[0][1]*intersectionArray[1][3]*intersectionArray[3][0]+intersectionArray[0][3]*intersectionArray[1][0]*intersectionArray[3][1]-intersectionArray[0][0]*intersectionArray[1][3]*intersectionArray[3][1]-intersectionArray[0][1]*intersectionArray[1][0]*intersectionArray[3][3]-intersectionArray[0][3]*intersectionArray[1][1]*intersectionArray[3][0];
            result.intersectionArray[2][3]= intersectionArray[0][0]*intersectionArray[1][3]*intersectionArray[2][1]+intersectionArray[0][1]*intersectionArray[1][0]*intersectionArray[2][3]+intersectionArray[0][3]*intersectionArray[1][1]*intersectionArray[2][3]-intersectionArray[0][0]*intersectionArray[1][1]*intersectionArray[2][3]-intersectionArray[0][1]*intersectionArray[1][3]*intersectionArray[2][0]-intersectionArray[0][3]*intersectionArray[1][0]*intersectionArray[2][1];
            
            result.intersectionArray[0][0] = intersectionArray[1][0]*intersectionArray[2][2]*intersectionArray[3][1]+intersectionArray[1][1]*intersectionArray[2][0]*intersectionArray[3][2]+intersectionArray[1][2]*intersectionArray[2][1]*intersectionArray[3][0]-intersectionArray[1][0]*intersectionArray[2][1]*intersectionArray[3][2]-intersectionArray[1][1]*intersectionArray[2][2]*intersectionArray[3][0]-intersectionArray[1][2]*intersectionArray[2][0]*intersectionArray[3][1];
            result.intersectionArray[0][1] = intersectionArray[0][0]*intersectionArray[2][1]*intersectionArray[3][2]+intersectionArray[0][1]*intersectionArray[2][2]*intersectionArray[3][0]+intersectionArray[0][2]*intersectionArray[2][0]*intersectionArray[3][1]-intersectionArray[0][0]*intersectionArray[2][2]*intersectionArray[3][1]-intersectionArray[0][1]*intersectionArray[2][0]*intersectionArray[3][2]-intersectionArray[0][2]*intersectionArray[2][1]*intersectionArray[3][0];
            result.intersectionArray[0][2]= intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[3][1]+intersectionArray[0][1]*intersectionArray[1][0]*intersectionArray[3][2]+intersectionArray[0][2]*intersectionArray[1][1]*intersectionArray[3][0]-intersectionArray[0][0]*intersectionArray[1][1]*intersectionArray[3][2]-intersectionArray[0][1]*intersectionArray[1][2]*intersectionArray[3][0]-intersectionArray[0][2]*intersectionArray[1][0]*intersectionArray[3][1];
            result.intersectionArray[0][3]= intersectionArray[0][0]*intersectionArray[1][1]*intersectionArray[2][2]+intersectionArray[0][1]*intersectionArray[1][2]*intersectionArray[2][0]+intersectionArray[0][2]*intersectionArray[1][0]*intersectionArray[2][1]-intersectionArray[0][0]*intersectionArray[1][2]*intersectionArray[2][1]-intersectionArray[0][1]*intersectionArray[1][0]*intersectionArray[2][2]-intersectionArray[0][2]*intersectionArray[1][1]*intersectionArray[2][0];
            //return result;


            float d = this->determinant();
        
            d = 1/d;
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    result.intersectionArray[i][j] *=d ;
                }
            }
           
            return result;
        }
       
        
        void operator = (Matrix copy)
        {
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    this->intersectionArray[i][j]=copy.intersectionArray[i][j];
                }
            }
        }
        

         void print()
        {   
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                     std::cout << intersectionArray[i][j]<< std::endl;
                } 
            }
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

        //face constructors
        Face(): v0_id(0),  v1_id(0), v2_id(0) {}
        Face(int a, int b, int c): v0_id(a),  v1_id(b), v2_id(c) {}

        void operator = (Face f) { v0_id=f.v0_id;  v1_id=f.v1_id; v2_id=f.v2_id; }
    };

struct Triangle
    {
        int material_id;
        Face indices;

        // the trianles helpers
        Vec3f v0, v1, v2 ;
        Vec3f n;
        Vec3f unitNormal;

        Material material;

        Triangle(): material_id(0) , indices(Face(0,0,0)) {}
        Triangle(int m, Face i): material_id(m) , indices(i) {}

        void computeNormal()
        {
            n = (v1-v0).cross(v2-v0);
            unitNormal = n.normalize();
        }
 
        float computeDet(float a,float b,float c,float d,float e,float f,float g,float h,float i)
        {
            float det = a*(e*i-h*f) + b*(g*f-d*i) + c*(d*h-e*g);
            return  det ;
        }

        bool isIntersect(Ray& ray, float& t)
        {
            // check if point is in triangle's plane
            if(unitNormal.dot(ray.d) == 0) 
            {
                return false;
            }

            float detA    = computeDet(v0.x-v1.x, v0.x-v2.x, ray.d.x,
                                         v0.y-v1.y, v0.y-v2.y, ray.d.y,
                                         v0.z-v1.z, v0.z-v2.z, ray.d.z);

            float detAlpha = computeDet(v0.x-ray.o.x, v0.x-v2.x, ray.d.x,
                                         v0.y-ray.o.y, v0.y-v2.y, ray.d.y, 
                                         v0.z-ray.o.z, v0.z-v2.z, ray.d.z);

            float detBeta = computeDet(v0.x-v1.x, v0.x-ray.o.x, ray.d.x, 
                                        v0.y-v1.y, v0.y-ray.o.y, ray.d.y,
                                        v0.z-v1.z, v0.z-ray.o.z, ray.d.z);

            float detT     = computeDet(v0.x-v1.x, v0.x-v2.x, v0.x-ray.o.x,
                                         v0.y-v1.y, v0.y-v2.y, v0.y-ray.o.y,
                                         v0.z-v1.z, v0.z-v2.z, v0.z-ray.o.z);

            float alpha = detAlpha / detA;
            float beta  = detBeta / detA;
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
        std::vector<Triangle> triangles;

    };

    

    struct Sphere
    {
        int material_id;
        int center_vertex_id;
        float radius;
        
        Vec3f sCenter;
        Vec3f unitNormal;
        Material sMaterial;
        
        //to compute normal line
        void computeNormal(Vec3f v)
        {
            unitNormal = (v-sCenter)/radius;
        }   

        //check intersction with the ray
        bool isIntersect(Ray& ray, float& t)
        {
           
            Vec3f c = this->sCenter;
            Vec3f o = ray.o;
            Vec3f d = ray.d;
            Vec3f oc = o-c;

            float R = this->radius;
            float B = 2*(d.dot(oc));
            float C = oc.dot(oc) - R*R ;
            float D = B*B - 4*C;

            if(D < -0.00001f)
            {   return false;}

            else
            {               
                float t1 = (-B - sqrt(D)) / 2;
                float t2 = (-B + sqrt(D)) / 2;
                
                if(t1<t2)
                {
                    t = t1;
                }
                else
                {
                    t = t2;
                }
                
                if(t<=0)
                {
                    return 0; 
                } 
               
                Vec3f rt = o + d*t ; // intersection formula
                computeNormal(rt);
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
<<<<<<< HEAD
		bool isIntersected(Ray ray, float& t, Material& imat, Vec3f& un);
		Vec3i computeAmbientLight(Ray ray, float& t, Material& material, Vec3f& un);
=======
        Vec3i computeShadow(Ray ray, float t, Vec3f n, Material material,int maxRec);
		bool isIntersected(Ray ray, float& t, Material& imat, Vec3f& un);
		

        
>>>>>>> 7bafb10aa13de0f7c8b9cf833df1e33240e7eace

    };

}

#endif
