//KITAA
#include "parser.h"
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include <iostream>
#define INFINTY numeric_limits<float>::infinity();

using namespace std;
using namespace parser;

void parser::Scene::loadFromXml(const std::string& filepath)
{
    tinyxml2::XMLDocument file;
    std::stringstream stream;

    auto res = file.LoadFile(filepath.c_str());
    if (res)
    {
        throw std::runtime_error("Error: The xml file cannot be loaded.");
    }

    auto root = file.FirstChild();
    if (!root)
    {
        throw std::runtime_error("Error: Root is not found.");
    }

    //Get BackgroundColor
    auto element = root->FirstChildElement("BackgroundColor");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0 0 0" << std::endl;
    }
    stream >> background_color.x >> background_color.y >> background_color.z;

    //Get ShadowRayEpsilon
    element = root->FirstChildElement("ShadowRayEpsilon");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0.001" << std::endl;
    }
    stream >> shadow_ray_epsilon;

    //Get MaxRecursionDepth
    element = root->FirstChildElement("MaxRecursionDepth");
    if (element)
    {
        stream << element->GetText() << std::endl;
    }
    else
    {
        stream << "0" << std::endl;
    }
    stream >> max_recursion_depth;

    //Get Cameras
    element = root->FirstChildElement("Cameras");
    element = element->FirstChildElement("Camera");
    Camera camera;
    while (element)
    {
        auto child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Gaze");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Up");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearPlane");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("NearDistance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageResolution");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("ImageName");
        stream << child->GetText() << std::endl;

        stream >> camera.position.x >> camera.position.y >> camera.position.z;
        stream >> camera.gaze.x >> camera.gaze.y >> camera.gaze.z;
        stream >> camera.up.x >> camera.up.y >> camera.up.z;
        stream >> camera.near_plane.x >> camera.near_plane.y >> camera.near_plane.z >> camera.near_plane.w;
        stream >> camera.near_distance;
        stream >> camera.image_width >> camera.image_height;
        stream >> camera.image_name;

        cameras.push_back(camera);
        element = element->NextSiblingElement("Camera");
    }

    //Get Lights
    element = root->FirstChildElement("Lights");
    auto child = element->FirstChildElement("AmbientLight");
    stream << child->GetText() << std::endl;
    stream >> ambient_light.x >> ambient_light.y >> ambient_light.z;
    element = element->FirstChildElement("PointLight");
    PointLight point_light;
    while (element)
    {
        child = element->FirstChildElement("Position");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("Intensity");
        stream << child->GetText() << std::endl;

        stream >> point_light.position.x >> point_light.position.y >> point_light.position.z;
        stream >> point_light.intensity.x >> point_light.intensity.y >> point_light.intensity.z;

        point_lights.push_back(point_light);
        element = element->NextSiblingElement("PointLight");
    }

    //Get Materials
    element = root->FirstChildElement("Materials");
    element = element->FirstChildElement("Material");
    Material material;
    while (element)
    {
        child = element->FirstChildElement("AmbientReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("DiffuseReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("SpecularReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("MirrorReflectance");
        stream << child->GetText() << std::endl;
        child = element->FirstChildElement("PhongExponent");
        stream << child->GetText() << std::endl;

        stream >> material.ambient.x >> material.ambient.y >> material.ambient.z;
        stream >> material.diffuse.x >> material.diffuse.y >> material.diffuse.z;
        stream >> material.specular.x >> material.specular.y >> material.specular.z;
        stream >> material.mirror.x >> material.mirror.y >> material.mirror.z;
        stream >> material.phong_exponent;

        materials.push_back(material);
        element = element->NextSiblingElement("Material");
    }

    //Get VertexData
    element = root->FirstChildElement("VertexData");
    stream << element->GetText() << std::endl;
    Vec3f vertex;
    while (!(stream >> vertex.x).eof())
    {
        stream >> vertex.y >> vertex.z;
        vertex_data.push_back(vertex);
    }
    stream.clear();

    //Get Meshes
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Mesh");
    Mesh mesh;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> mesh.material_id;

        child = element->FirstChildElement("Faces");
        stream << child->GetText() << std::endl;
        Face face;
        while (!(stream >> face.v0_id).eof())
        {
            stream >> face.v1_id >> face.v2_id;
            mesh.faces.push_back(face);
        }
        stream.clear();

        meshes.push_back(mesh);
        mesh.faces.clear();
        element = element->NextSiblingElement("Mesh");
    }
    stream.clear();

    //Get Triangles
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Triangle");
    Triangle triangle;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> triangle.material_id;

        child = element->FirstChildElement("Indices");
        stream << child->GetText() << std::endl;
        stream >> triangle.indices.v0_id >> triangle.indices.v1_id >> triangle.indices.v2_id;

        triangles.push_back(triangle);
        element = element->NextSiblingElement("Triangle");
    }

    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("Sphere");
    Sphere sphere;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphere.radius;

        spheres.push_back(sphere);
        element = element->NextSiblingElement("Sphere");
    }


}

void parser::Scene::init()
{

    point_lights_size = point_lights.size();

    // init spheres
    spheres_size = spheres.size() ; 
    for(int sph=0 ; sph < spheres_size ; sph++)
    {
        //spheres[sph].compute_transform(scalingMatrices, translationMatrices, rotationMatrices);
        spheres[sph].center = vertex_data[spheres[sph].center_vertex_id-1];
        spheres[sph].mat = materials[spheres[sph].material_id-1] ;
    }

    // init triangles
    triangles_size = triangles.size() ;
    for(int tri=0 ; tri < triangles_size; tri++)
    {               
        triangles[tri].vertex0 = vertex_data[triangles[tri].indices.v0_id-1];
        triangles[tri].vertex1 = vertex_data[triangles[tri].indices.v1_id-1];
        triangles[tri].vertex2 = vertex_data[triangles[tri].indices.v2_id-1];
        triangles[tri].compute_normal();
        //triangles[tri].transform_triangle(scalingMatrices, translationMatrices, rotationMatrices); //hw2

        triangles[tri].mat = materials[triangles[tri].material_id-1];
    }

    // init meshes
    meshes_size = meshes.size() ;
    for(int msh=0; msh < meshes_size; msh++)
    {
        Material msh_mat =  materials[meshes[msh].material_id-1];

        for(int f=0; f < meshes[msh].faces.size(); f++)
        {
            Triangle tri(meshes[msh].material_id, meshes[msh].faces[f]);

            tri.vertex0 = vertex_data[tri.indices.v0_id-1];
            tri.vertex1 = vertex_data[tri.indices.v1_id-1];
            tri.vertex2 = vertex_data[tri.indices.v2_id-1]; 
            tri.compute_normal();
            tri.mat = msh_mat;

            meshes[msh].mtriangles.push_back(tri);  
        }

     }

   
}   



//check if ray-spehere intersects or triangle-ray intersects
bool parser::Scene::isIntersected(Ray ray, float& t, Material& material, Vec3f& unitNormalVector)
{
	float tMin = INFINITY;
	int sphereNumber = spheres.size();
	int triangleNumber = triangles.size();
	int meshNumber = meshes.size();
	int meshInstanceNumber = meshes.size();

	for (int sphereIndex = 0; sphereIndex < sphereNumber; sphereIndex++) //size from initScene
	{
		if (spheres[sphereIndex].is_intersect(ray, t) && t<tMin)
		{
			tMin = t;
			material = spheres[sphereIndex].mat;
			//sphereIndexeres[sphereIndex].compute_normal(ray.origin+ray.direction*t);
			unitNormalVector = spheres[sphereIndex].unit_normal;
		}
	}

	for (int triangleIndex = 0; triangleIndex < triangleNumber; triangleIndex++)
	{
		if (triangles[triangleIndex].is_intersect(ray, t) && t<tMin)
		{
			tMin = t;
			material = triangles[triangleIndex].mat;
			unitNormalVector = triangles[triangleIndex].unit_normal;
		}
	}

	for (int meshIndex = 0; meshIndex < meshNumber; meshIndex++)
	{
		for (int f = 0; f < meshes[meshIndex].faces.size(); f++)
		{
			if (meshes[meshIndex].mtriangles[f].is_intersect(ray, t) && t<tMin)
			{
				tMin = t;
				material = meshes[meshIndex].mtriangles[f].mat;
				unitNormalVector = meshes[meshIndex].mtriangles[f].unit_normal;
			}
		}
	}

	t = tMin;

	if (tMin != INFINITY)
		return true;
	else
		return false;
}


float clamping(float intensityValue)
{
    if (intensityValue > 255.0)
        return 255.0;
    else
        return intensityValue;
}


 //IF USING THESE HERE,ADD THEIR DECLARATION IN PARSER HEADER FILE!!!
Vec3i parser::Scene::computeAmbientLight(Ray ray, float& t, Material& material, Vec3f& unitNormal, int& count)
{
	Vec3f colorAmbient;
	colorAmbient.x = (int)clamping(ambient_light.x * material.ambient.x);
	colorAmbient.y = (int)clamping(ambient_light.y * material.ambient.y);
	colorAmbient.z = (int)clamping(ambient_light.z * material.ambient.z);

}

//*****************************CALCULATE RAY-TRIANGLE INTERSECTION AND RAY-SPHERE INTERSECTION***********************************************////////////USE BLITZ LIBRARY FOR MATRICE OPERATIONS