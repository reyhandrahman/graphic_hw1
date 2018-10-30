#include "parser.h"
#include "tinyxml2.h"
#include <sstream>
#include <stdexcept>
#include <iostream>

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

    //Get sphereIndexeres
    element = root->FirstChildElement("Objects");
    element = element->FirstChildElement("sphereIndexere");
    sphereIndexere sphereIndexere;
    while (element)
    {
        child = element->FirstChildElement("Material");
        stream << child->GetText() << std::endl;
        stream >> sphereIndexere.material_id;

        child = element->FirstChildElement("Center");
        stream << child->GetText() << std::endl;
        stream >> sphereIndexere.center_vertex_id;

        child = element->FirstChildElement("Radius");
        stream << child->GetText() << std::endl;
        stream >> sphereIndexere.radius;

        sphereIndexeres.push_back(sphereIndexere);
        element = element->NextSiblingElement("sphereIndexere");
    }
}
bool parser::Scene::isIntersected(Ray ray, float& t, Material& imat, Vec3f& un)
{
	float tMin = numeric<float>::infinity();
	int sphereNumber = spheres.size();
	int triangleNumber = triangles.size();
	int meshNumber = meshes.size();
	int meshInstanceNumber = meshIstances.size();
	for (int sphereIndex = 0; sphereIndex < sphereNumber; sphereIndex++) //soze from initScene
	{
		if (sphereIndexeres[sphereIndex].is_intersect(ray, t) && t<tMin)
		{
			tMin = t;
			imat = sphereIndexeres[sphereIndex].mat;
			//sphereIndexeres[sphereIndex].compute_normal(ray.origin+ray.direction*t);
			un = sphereIndexeres[sphereIndex].unit_normal;
		}
	}

	for (int triangleIndex = 0; triangleIndex < triangleNumber; triangleIndex++)
	{
		if (triangles[triangleIndex].is_intersect(ray, t) && t<tMin)
		{
			tMin = t;
			imat = triangles[triangleIndex].mat;
			un = triangles[triangleIndex].unit_normal;
		}
	}

	for (int meshIndex = 0; meshIndex < meshes_size; meshIndex++)
	{
		for (int f = 0; f < meshes[meshIndex].faces.size(); f++)
		{
			if (meshes[meshIndex].mtriangles[f].is_intersect(ray, t) && t<tMin)
			{
				tMin = t;
				imat = meshes[meshIndex].mtriangles[f].mat;
				un = meshes[meshIndex].mtriangles[f].unit_normal;
			}
		}
	}

	for (int meshInstanceIndex = 0; meshInstanceIndex < meshInstance_size; meshInstanceIndex++)
	{
		for (int f = 0; f < meshInstances[meshInstanceIndex].baseMesh.faces.size(); f++)
		{
			if (meshInstances[meshInstanceIndex].baseMesh.mtriangles[f].is_intersect(ray, t) && t<tMin)
			{
				tMin = t;
				imat = meshInstances[meshInstanceIndex].baseMesh.mtriangles[f].mat;
				un = meshInstances[meshInstanceIndex].baseMesh.mtriangles[f].unit_normal;
			}
		}
	}

	t = tMin;

	if (tMin != numeric_limits<float>::infinity())
		return true;
	else
		return false;
}

 //IF USING THESE HERE,ADD THEIR DECLARATION IN PARSER HEADER FILE!!!
Vec3i parser::Scene::computeAmbientLight(Ray ray, float& t, Material& material, Vec3f& un)
{
	Vec3f colorAmbient;
	colorAmbient.x = (int)clamping(ambient_light.x * material.ambient.x);
	colorAmbient.y = (int)clamping(ambient_light.y * material.ambient.y);
	colorAmbient.z = (int)clamping(ambient_light.z * material.ambient.z);

}

float clamping(float intensityValue)
{
	if (intensityValue > 255.0)
		return 255.0;
	else
		return intensityValue;
}

//*****************************CALCULATE RAY-TRIANGLE INTERSECTION AND RAY-SPHERE INTERSECTION***********************************************////////////USE BLITZ LIBRARY FOR MATRICE OPERATIONS