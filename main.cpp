// ini punya kitaaa
#include <iostream>
#include "parser.h"
#include "ppm.h"

using namespace parser;
using namespace std;

typedef unsigned char RGB[3];

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    // The code below creates a test pattern and writes
    // it to a PPM file to demonstrate the usage of the
    // ppm_write function.

	

	//prepare everything, to ease and fasten the computing
	scene.init();
	
	//for each camera
	for (int cameraIndex = 0; cameraIndex < scene.cameras.size(); cameraIndex++)
	{
		int imageWidth = scene.cameras[cameraIndex].image_width;
		int imageHeight = scene.cameras[cameraIndex].image_height;

		Vec3f e = scene.cameras[cameraIndex].position;
		Vec3f w = scene.cameras[cameraIndex].gaze * (-1);
		Vec3f v = scene.cameras[cameraIndex].up;
		Vec3f u = v.cross(w); //cross product

		w = w.normalize();
		v = v.normalize();
		u = u.normalize();

		float distance = scene.cameras[cameraIndex].near_distance;

		float b = scene.cameras[cameraIndex].near_plane.z; 
		float l = scene.cameras[cameraIndex].near_plane.x;
		float t = scene.cameras[cameraIndex].near_plane.w;
		float r = scene.cameras[cameraIndex].near_plane.y;

		Vec3f s;
		float su, sv;

		Vec3f m = e + ((w*(-1))*distance); //REPLACE W*-1 WITH GAZE 
		Vec3f q = m + u * l + v * t;

		

		float rightToLeftUnit = (r - l) / imageWidth;
		float topToBottomUnit = (t - b) / imageHeight;

		unsigned char* image = new unsigned char[imageWidth * imageHeight * 3];

		int index = 0;


		//for each pixel
		for (int j = 0; j < imageHeight; ++j)
		{
			for (int i = 0; i < imageWidth; ++i)
			{
				//compute s
				su = rightToLeftUnit * (i + 0.50);
				sv = topToBottomUnit * (j + 0.50);
				s = q + (u * su) - (v * sv);

				Vec3f rayDirection = (s - e).normalize(); 
				Ray ray(e, rayDirection); 

				//new variables
				float t; 
				Material material;
				Vec3f un;

				//check ray-sphere intersection or ray-triangle intersection
				//load the materials here

				if (scene.isIntersected(ray, t, material, un)) 
				{
					//Vec3i color = scene.computeAmbientLight(ray, t, material, un, scene.max_recursion_depth);
					Vec3i color    = scene.computeShadow(ray, t, un, material, scene.max_recursion_depth);
					image[index++] = color.x;
					image[index++] = color.y;
					image[index++] = color.z;
				}
				else
				{
					image[index++] = scene.background_color.x;
					image[index++] = scene.background_color.y;
					image[index++] = scene.background_color.z;
				}
			}
		}
		//*****************************************************************
    

  
	
    	write_ppm(argv[2], image, imageWidth, imageHeight);
    }

}
