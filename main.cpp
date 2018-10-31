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

	//**************************************************
	
	for (int cameraIndex = 0; cameraIndex < scene.cameras.size(); cameraIndex++)
	{
		int imageWidth = scene.cameras[cameraIndex].image_width;
		int imageHeight = scene.cameras[cameraIndex].image_height;

		Vec3f e = scene.cameras[cameraIndex].position;
		Vec3f w = scene.cameras[cameraIndex].gaze * (-1);
		Vec3f v = scene.cameras[cameraIndex].up;
		Vec3f u = v.cross(w);

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


		int width = 640, height = 480;
		//for each pixel
		for (int j = 0; j < height; ++j)
		{
			for (int i = 0; i < width; ++i)
			{
				//compute s
				su = rightToLeftUnit * (i + 0.50);
				sv = topToBottomUnit * (j + 0.50);
				s = q + (u * su) - (v * sv);

				Vec3f rayDirection = (s - e).normalize(); 
				//se = se.normalize();

				Ray ray(e, rayDirection); //DEFINE WHAT'S A RAY

				float t;
				Material material;
				Vec3f un;

				// give color value
				if (scene.isIntersected(ray, t, material, un)) 
				{
					Vec3i color = scene.computeAmbientLight(ray, t, material, un, scene.max_recursion_depth); //calculate_color define
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
    

    // const RGB BAR_COLOR[8] =
    // {
    //     { 255, 255, 255 },  // 100% White
    //     { 255, 255,   0 },  // Yellow
    //     {   0, 255, 255 },  // Cyan
    //     {   0, 255,   0 },  // Green
    //     { 255,   0, 255 },  // Magenta
    //     { 255,   0,   0 },  // Red
    //     {   0,   0, 255 },  // Blue
    //     {   0,   0,   0 },  // Black
    // };

    // int width = 640, height = 480;
    // int columnWidth = width / 8;

    // unsigned char* image = new unsigned char [width * height * 3];

    // int i = 0;
    // for (int y = 0; y < height; ++y)
    // {
    //     for (int x = 0; x < width; ++x)
    //     {
    //         int colIdx = x / columnWidth;
    //         image[i++] = scene.background_color.x;
    //         image[i++] = scene.background_color.y;
    //         image[i++] = scene.background_color.z;
    //     }
    // }

    // //try fot triangle
    // for (int y = 0; y < height; ++y)
    // {
    //     for (int x = 0; x < width; ++x)
    //     {
    //         int colIdx = x / columnWidth;
    //         image[i++] = scene.background_color.x;
    //         image[i++] = scene.background_color.y;
    //         image[i++] = scene.background_color.z;
    //     }
    // }

	
    	write_ppm(argv[2], image, width, height);
    }

}
/* //PASS A SCENE AS AN ARGUMENT OR HAVE IT IN PARSER-time inefficient???
Vec3i computeAmbientLight(Ray ray, float& t, Material& material, Vec3f& un)
{
	Vec3f colorAmbient;
	colorAmbient.x = (int)ambient_light.x * material.ambient.x;
	colorAmbient.y = (int)ambient_light.y * material.ambient.y;
	colorAmbient.z = (int)ambient_light.z * material.ambient.z;

}

float clamping(float intensityValue)
{
	if (intensityValue > 255.0)
		return 255.0;
	else
		return intensityValue;
}
*/