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

    //init sphere
	for(int s=0 ; s < scene.spheres.size() ; s++)
    {
        //id minus 1 so it starts from 0
        scene.spheres[s].sCenter = scene.vertex_data[scene.spheres[s].center_vertex_id-1];
        scene.spheres[s].sMaterial = scene.materials[scene.spheres[s].material_id-1] ;
    }

    // init triangles
    for(int t=0 ; t < scene.triangles.size(); t++)
    {               
        scene.triangles[t].v0 = scene.vertex_data[scene.triangles[t].indices.v0_id-1];
        scene.triangles[t].v1 = scene.vertex_data[scene.triangles[t].indices.v1_id-1];
        scene.triangles[t].v2 = scene.vertex_data[scene.triangles[t].indices.v2_id-1];
        scene.triangles[t].computeNormal();

        scene.triangles[t].material = scene.materials[scene.triangles[t].material_id-1];
    }

    // init meshes
    for(int msh=0; msh < scene.meshes.size(); msh++)
    {
        Material msh_mat =  scene.materials[scene.meshes[msh].material_id-1];

        for(int f=0; f < scene.meshes[msh].faces.size(); f++)
        {
            Triangle tri(scene.meshes[msh].material_id, scene.meshes[msh].faces[f]);

            tri.v0 = scene.vertex_data[tri.indices.v0_id-1];
            tri.v1 = scene.vertex_data[tri.indices.v1_id-1];
            tri.v2 = scene.vertex_data[tri.indices.v2_id-1]; 
            tri.computeNormal();
            tri.material = msh_mat;

            scene.meshes[msh].triangles.push_back(tri);  
        }

     }
	//for each camera
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

		int i = 0;
		for (int y = 0; y < imageHeight; ++y)
		{
			for (int x = 0; x < imageWidth; ++x)
			{
				//compute s
				su = rightToLeftUnit * (x + 0.50);
				sv = topToBottomUnit * (y + 0.50);
				s = q + (u * su) - (v * sv);

				Vec3f rayDirection = (s - e).normalize(); 
				Ray ray(e, rayDirection); 

				float t;
				Material material;
				Vec3f un;

				//check ray-sphere intersection or ray-triangle intersection
				if (scene.isIntersected(ray, t, material, un)) 
				{
					//Vec3i color = scene.computeAmbientLight(ray, t, material, un, scene.max_recursion_depth);
					Vec3i color  = scene.computeShadow(ray, t, un, material, scene.max_recursion_depth);
					image[i++] = color.x;
					image[i++] = color.y;
					image[i++] = color.z;
				}
				else
				{
					image[i++] = scene.background_color.x;
					image[i++] = scene.background_color.y;
					image[i++] = scene.background_color.z;
				}
			}
		}
		
	
    	write_ppm(argv[2], image, imageWidth, imageHeight);
    }

}
