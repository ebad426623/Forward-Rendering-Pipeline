#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>

#include "tinyxml2.h"
#include "Triangle.h"
#include "Helpers.h"
#include "Scene.h"

using namespace tinyxml2;
using namespace std;



std::vector<std::vector<double>> depth;

/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *xmlElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *rootNode = xmlDoc.FirstChild();

	// read background color
	xmlElement = rootNode->FirstChildElement("BackgroundColor");
	str = xmlElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	xmlElement = rootNode->FirstChildElement("Culling");
	if (xmlElement != NULL)
	{
		str = xmlElement->GetText();

		if (strcmp(str, "enabled") == 0)
		{
			this->cullingEnabled = true;
		}
		else
		{
			this->cullingEnabled = false;
		}
	}

	// read cameras
	xmlElement = rootNode->FirstChildElement("Cameras");
	XMLElement *camElement = xmlElement->FirstChildElement("Camera");
	XMLElement *camFieldElement;
	while (camElement != NULL)
	{
		Camera *camera = new Camera();

		camElement->QueryIntAttribute("id", &camera->cameraId);

		// read projection type
		str = camElement->Attribute("type");

		if (strcmp(str, "orthographic") == 0)
		{
			camera->projectionType = ORTOGRAPHIC_PROJECTION;
		}
		else
		{
			camera->projectionType = PERSPECTIVE_PROJECTION;
		}

		camFieldElement = camElement->FirstChildElement("Position");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->position.x, &camera->position.y, &camera->position.z);

		camFieldElement = camElement->FirstChildElement("Gaze");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->gaze.x, &camera->gaze.y, &camera->gaze.z);

		camFieldElement = camElement->FirstChildElement("Up");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf", &camera->v.x, &camera->v.y, &camera->v.z);

		camera->gaze = normalizeVec3(camera->gaze);
		camera->u = crossProductVec3(camera->gaze, camera->v);
		camera->u = normalizeVec3(camera->u);

		camera->w = inverseVec3(camera->gaze);
		camera->v = crossProductVec3(camera->u, camera->gaze);
		camera->v = normalizeVec3(camera->v);

		camFieldElement = camElement->FirstChildElement("ImagePlane");
		str = camFieldElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			   &camera->left, &camera->right, &camera->bottom, &camera->top,
			   &camera->near, &camera->far, &camera->horRes, &camera->verRes);

		camFieldElement = camElement->FirstChildElement("OutputName");
		str = camFieldElement->GetText();
		camera->outputFilename = string(str);

		this->cameras.push_back(camera);

		camElement = camElement->NextSiblingElement("Camera");
	}

	// read vertices
	xmlElement = rootNode->FirstChildElement("Vertices");
	XMLElement *vertexElement = xmlElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (vertexElement != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

		vertex->colorId = vertexId;

		str = vertexElement->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = vertexElement->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		this->vertices.push_back(vertex);
		this->colorsOfVertices.push_back(color);

		vertexElement = vertexElement->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	xmlElement = rootNode->FirstChildElement("Translations");
	XMLElement *translationElement = xmlElement->FirstChildElement("Translation");
	while (translationElement != NULL)
	{
		Translation *translation = new Translation();

		translationElement->QueryIntAttribute("id", &translation->translationId);

		str = translationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		this->translations.push_back(translation);

		translationElement = translationElement->NextSiblingElement("Translation");
	}

	// read scalings
	xmlElement = rootNode->FirstChildElement("Scalings");
	XMLElement *scalingElement = xmlElement->FirstChildElement("Scaling");
	while (scalingElement != NULL)
	{
		Scaling *scaling = new Scaling();

		scalingElement->QueryIntAttribute("id", &scaling->scalingId);
		str = scalingElement->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		this->scalings.push_back(scaling);

		scalingElement = scalingElement->NextSiblingElement("Scaling");
	}

	// read rotations
	xmlElement = rootNode->FirstChildElement("Rotations");
	XMLElement *rotationElement = xmlElement->FirstChildElement("Rotation");
	while (rotationElement != NULL)
	{
		Rotation *rotation = new Rotation();

		rotationElement->QueryIntAttribute("id", &rotation->rotationId);
		str = rotationElement->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		this->rotations.push_back(rotation);

		rotationElement = rotationElement->NextSiblingElement("Rotation");
	}

	// read meshes
	xmlElement = rootNode->FirstChildElement("Meshes");

	XMLElement *meshElement = xmlElement->FirstChildElement("Mesh");
	while (meshElement != NULL)
	{
		Mesh *mesh = new Mesh();

		meshElement->QueryIntAttribute("id", &mesh->meshId);

		// read projection type
		str = meshElement->Attribute("type");

		if (strcmp(str, "wireframe") == 0)
		{
			mesh->type = WIREFRAME_MESH;
		}
		else
		{
			mesh->type = SOLID_MESH;
		}

		// read mesh transformations
		XMLElement *meshTransformationsElement = meshElement->FirstChildElement("Transformations");
		XMLElement *meshTransformationElement = meshTransformationsElement->FirstChildElement("Transformation");

		while (meshTransformationElement != NULL)
		{
			char transformationType;
			int transformationId;

			str = meshTransformationElement->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			mesh->transformationTypes.push_back(transformationType);
			mesh->transformationIds.push_back(transformationId);

			meshTransformationElement = meshTransformationElement->NextSiblingElement("Transformation");
		}

		mesh->numberOfTransformations = mesh->transformationIds.size();

		// read mesh faces
		char *row;
		char *cloneStr;
		int v1, v2, v3;
		XMLElement *meshFacesElement = meshElement->FirstChildElement("Faces");
		str = meshFacesElement->GetText();
		cloneStr = strdup(str);

		row = strtok(cloneStr, "\n");
		while (row != NULL)
		{
			int result = sscanf(row, "%d %d %d", &v1, &v2, &v3);

			if (result != EOF)
			{
				mesh->triangles.push_back(Triangle(v1, v2, v3));
			}
			row = strtok(NULL, "\n");
		}
		mesh->numberOfTriangles = mesh->triangles.size();
		this->meshes.push_back(mesh);

		meshElement = meshElement->NextSiblingElement("Mesh");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera *camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;
			vector<double> rowOfDepths;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
				rowOfDepths.push_back(1.01);
			}

			this->image.push_back(rowOfColors);
			depth.push_back(rowOfDepths);
		}
	}
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
				depth[i][j] = 1.01;
				depth[i][j] = 1.01;
				depth[i][j] = 1.01;
			}
		}
	}
}

/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}


/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera *camera)
{
	ofstream fout;

	fout.open(camera->outputFilename.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFilename << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].g) << " "
				 << makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "./magick " + ppmFileName + " " + ppmFileName + ".png";
		int a = system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick " + ppmFileName + " " + ppmFileName + ".png";
		int a = system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}


/*
	Transformations, clipping, culling, rasterization are done here.
*/


struct Vec3_helper
{
	int x, y;
	double z, r, g, b;
};

struct Vec4_helper
{
	double x, y, z, t, r, g, b;
};

int findMin(int a, int b, int c) {
    if (a < b) {
        return (a < c) ? a : c;
    } else {
        return (b < c) ? b : c;
    }
}

int findMax(int a, int b, int c) {
    if (a > b) {
        return (a > c) ? a : c;
    } else {
        return (b > c) ? b : c;
    }
}

void camera_projection_transformation(Camera &camera, Matrix4 &projection_trans) {

	projection_trans.values[0][3] = inverseVec3(camera.position).x;
	projection_trans.values[1][3] = inverseVec3(camera.position).y;
	projection_trans.values[2][3] = inverseVec3(camera.position).z;


    Matrix4 projection_trans_1 = getIdentityMatrix();

	projection_trans_1.values[0][0] = camera.u.x;
	projection_trans_1.values[0][1] = camera.u.y;
	projection_trans_1.values[0][2] = camera.u.z;

	projection_trans_1.values[1][0] = camera.v.x;
	projection_trans_1.values[1][1] = camera.v.y;
	projection_trans_1.values[1][2] = camera.v.z;

	projection_trans_1.values[2][0] = camera.w.x;
	projection_trans_1.values[2][1] = camera.w.y;
	projection_trans_1.values[2][2] = camera.w.z;


    projection_trans = multiplyMatrixWithMatrix(projection_trans_1, projection_trans);


    double right_left = camera.right - camera.left;
    double top_bottom = camera.top - camera.bottom;
    double far_near = camera.far - camera.near;


    projection_trans_1 = getIdentityMatrix();

    
	if (camera.projectionType) {
        projection_trans_1.values[0][0] = (2 * camera.near) / right_left;
        projection_trans_1.values[0][2] = (camera.right + camera.left) / right_left;
        projection_trans_1.values[1][1] = (2 * camera.near) / top_bottom;
        projection_trans_1.values[1][2] = (camera.top + camera.bottom) / top_bottom;
        projection_trans_1.values[2][2] = -(camera.far + camera.near) / far_near;
        projection_trans_1.values[2][3] = (-2 * camera.far * camera.near) / far_near;
        projection_trans_1.values[3][2] = -1;
        projection_trans_1.values[3][3] = 0;
    } else {
        projection_trans_1.values[0][0] = 2 / right_left;
        projection_trans_1.values[0][3] = -(camera.right + camera.left) / right_left;
        projection_trans_1.values[1][1] = 2 / top_bottom;
        projection_trans_1.values[1][3] = -(camera.top + camera.bottom) / top_bottom;
        projection_trans_1.values[2][2] = -2 / far_near;
        projection_trans_1.values[2][3] = -(camera.far + camera.near) / far_near;
    }

    projection_trans = multiplyMatrixWithMatrix(projection_trans_1, projection_trans);
}

void find_model_transformation(Scene &scene, Mesh &mesh, Matrix4 &model_trans) {
	
	int i = 0;
    
	while (i < mesh.numberOfTransformations) {
        Matrix4 model_trans_1 = getIdentityMatrix();

        if (mesh.transformationTypes[i] == 'r') { 
            Rotation rotation = *(scene.rotations[mesh.transformationIds[i] - 1]);
			Vec3 u, v, w;

			u.x = rotation.ux;
			u.y = rotation.uy;
			u.z = rotation.uz;

            if (rotation.ux <= rotation.uy && rotation.ux <= rotation.uz) {
				v.x = 0;
				v.y = -rotation.uz;
				v.z = rotation.uy;
			}
            else if (rotation.uy <= rotation.ux && rotation.uy <= rotation.uz) {
				v.x = -rotation.uz;
				v.y = 0;
				v.z = rotation.ux;
			}
            else {
				v.x = -rotation.uy;
				v.y = rotation.ux;
				v.z = 0;
			}



            normalizeVec3(u);
			model_trans_1.values[0][0] = u.x;
			model_trans_1.values[0][1] = u.y;
			model_trans_1.values[0][2] = u.z;

		
            normalizeVec3(v);
			model_trans_1.values[1][0] = v.x;
			model_trans_1.values[1][1] = v.y;
			model_trans_1.values[1][2] = v.z;



			w = crossProductVec3(u, v);
            normalizeVec3(w);
			model_trans_1.values[2][0] = w.x;
			model_trans_1.values[2][1] = w.y;
			model_trans_1.values[2][2] = w.z;


            model_trans = multiplyMatrixWithMatrix(model_trans_1, model_trans);
            model_trans_1 = getIdentityMatrix();
			double alpha = (rotation.angle * (3.14159265358979323846)) / 180;


			model_trans_1.values[1][1] = cos(alpha);
            model_trans_1.values[1][2] = -sin(alpha);
            model_trans_1.values[2][1] = sin(alpha);
            model_trans_1.values[2][2] = cos(alpha);


            model_trans = multiplyMatrixWithMatrix(model_trans_1, model_trans); 
            model_trans_1 = getIdentityMatrix();


			model_trans_1.values[0][0] = u.x;
			model_trans_1.values[1][0] = u.y;
			model_trans_1.values[2][0] = u.z;


			model_trans_1.values[0][1] = v.x;
			model_trans_1.values[1][1] = v.y;
			model_trans_1.values[2][1] = v.z;


			model_trans_1.values[0][2] = w.x;
			model_trans_1.values[1][2] = w.y;
			model_trans_1.values[2][2] = w.z;

            model_trans = multiplyMatrixWithMatrix(model_trans_1, model_trans);
		}

        else if (mesh.transformationTypes[i] == 't') { 
            Translation translation = *(scene.translations[mesh.transformationIds[i] - 1]);
            model_trans_1.values[0][3] = translation.tx;
            model_trans_1.values[1][3] = translation.ty;
            model_trans_1.values[2][3] = translation.tz;
            model_trans = multiplyMatrixWithMatrix(model_trans_1, model_trans);
        }

        else { 
            Scaling scale = *(scene.scalings[mesh.transformationIds[i] - 1]);
            model_trans_1.values[0][0] = scale.sx;
            model_trans_1.values[1][1] = scale.sy;
            model_trans_1.values[2][2] = scale.sz;
            model_trans = multiplyMatrixWithMatrix(model_trans_1, model_trans);
        }

        i++;
    }
}

void rasterize_mesh(Scene &scene, Camera &camera, vector<Vec4> &faces, Matrix4 &viewport_transformation) {
	vector<Vec3_helper> triangle;
	int i = 0;
	while(i < 3) {
		Vec4 temp1 = multiplyMatrixWithVec4(viewport_transformation, faces[i]);
		Vec3_helper temp2;
		Color c = *scene.colorsOfVertices[temp1.colorId - 1];
		temp2.x = static_cast<int>(temp1.x);
		temp2.y = static_cast<int>(temp1.y);
		temp2.z = temp1.z;
		temp2.r = c.r;
		temp2.g = c.g;
		temp2.b = c.b;
		triangle.push_back(temp2);
		i++;
	}
	
	Vec3_helper v[3];
	i = 0;

	while(i < 3) {
		v[i].x = static_cast<int>(triangle[i].x);
		v[i].y = static_cast<int>(triangle[i].y);
		v[i].z = triangle[i].z;
		v[i].r = triangle[i].r;
		v[i].g = triangle[i].g;
		v[i].b = triangle[i].b;
 		i++;
	}
	
	
	int xmin = findMin(v[0].x, v[1].x, v[2].x);
	int ymin = findMin(v[0].y, v[1].y, v[2].y);
	int xmax = findMax(v[0].x, v[1].x, v[2].x);
	int ymax = findMax(v[0].y, v[1].y, v[2].y);

	
	int x0 = v[0].x, y0 = v[0].y;
	int x1 = v[1].x, y1 = v[1].y;
	int x2 = v[2].x, y2 = v[2].y;

	double f01 = x2 * (y0 - y1) + y2 * (x1 - x0) + x0 * y1 - y0 * x1;
	double f12 = x0 * (y1 - y2) + y0 * (x2 - x1) + x1 * y2 - y1 * x2;
	double f20 = x1 * (y2 - y0) + y1 * (x0 - x2) + x2 * y0 - y2 * x0;


	int y = ymin;
	while (y <= ymax) {
		int x = xmin;
		while (x <= xmax) {

			if (f12 != 0) {
				int a1 = x * (v[1].y - v[2].y);
				int a2 = y * (v[2].x - v[1].x);
				int a3 = v[1].x * v[2].y - v[1].y * v[2].x;
				double alpha = static_cast<double>(a1 + a2 + a3) / f12;

				if (alpha >= 0) {
					if (f20 != 0) {
						int b1 = x * (v[2].y - v[0].y);
						int b2 = y * (v[0].x - v[2].x);
						int b3 = v[2].x * v[0].y - v[2].y * v[0].x;
						double beta = static_cast<double>(b1 + b2 + b3) / f20;

						if (beta >= 0) {
							if (f01 != 0) {
								int c1 = x * (v[0].y - v[1].y);
								int c2 = y * (v[1].x - v[0].x);
								int c3 = v[0].x * v[1].y - v[0].y * v[1].x;
								double gamma = static_cast<double>(c1 + c2 + c3) / f01;

								if (gamma >= 0) {
									double z = alpha * v[0].z + beta * v[1].z + gamma * v[2].z;

									if (x >= 0 && x < camera.horRes) {
										if (y >= 0 && y < camera.verRes) {
											if (z < depth[x][y]) {
												Color c;
												c.r = alpha * v[0].r + beta * v[1].r + gamma * v[2].r;
												c.g = alpha * v[0].g + beta * v[1].g + gamma * v[2].g;
												c.b = alpha * v[0].b + beta * v[1].b + gamma * v[2].b;

												scene.image[x][y].r = c.r;
												scene.image[x][y].g = c.g;
												scene.image[x][y].b = c.b;
												
												depth[x][y] = z;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			x++;
		}
		y++;
	}

}

bool is_visible(double d, double d1, double &enter, double &leave) {
    double t = d1 / d;

    if (d == 0 && d1 > 0) return false;

    if (d > 0) {
        if (t > leave) return false;

        enter = (t > enter) ? t : enter; 
    } 
	else if (d < 0) {
        if (t < enter) return false; 
		
        leave = (t < leave) ? t : leave; 
    }
	
	return true;
}

bool to_clip(Scene &scene, Camera &camera, Vec4 face0, Vec4 face1, vector<Vec4_helper> &temp) {
    double enter = 0, leave = 1;
    double delta_x = face1.x - face0.x;
	double delta_y = face1.y - face0.y;
	double delta_z = face1.z - face0.z;

    if (is_visible(delta_x, -1 - face0.x, enter, leave) &&
        is_visible(-delta_x, face0.x - 1, enter, leave) &&
        is_visible(delta_y, -1 - face0.y, enter, leave) &&
        is_visible(-delta_y, face0.y - 1, enter, leave) &&
        is_visible(delta_z, -1 - face0.z, enter, leave) &&
        is_visible(-delta_z, face0.z - 1, enter, leave)) {

        temp.clear();

        Color c0 = *scene.colorsOfVertices[face0.colorId - 1];
        Color c1 = *scene.colorsOfVertices[face1.colorId - 1];
        Vec4_helper new_face0, new_face1;
        
		new_face0 = (enter > 0) ? Vec4_helper{
            enter * delta_x + face0.x,
            enter * delta_y + face0.y,
            enter * delta_z + face0.z,
            face0.t,
            enter * c1.r + (1 - enter) * c0.r,
            enter * c1.g + (1 - enter) * c0.g,
            enter * c1.b + (1 - enter) * c0.b
        } : Vec4_helper{
			face0.x,
			face0.y,
			face0.z,
			face0.t,
			c0.r,
			c0.g,
			c0.b,
		};
        temp.push_back(new_face0);



        new_face1 = (leave < 1) ? Vec4_helper{
            leave * delta_x + face0.x,
            leave * delta_y + face0.y,
            leave * delta_z + face0.z,
            face1.t,
            leave * c1.r + (1 - leave) * c0.r ,
            leave * c1.g + (1 - leave) * c0.g ,
            leave * c1.b + (1 - leave) * c0.b
        } : Vec4_helper{
			face1.x,
			face1.y,
			face1.z,
			face1.t,
			c1.r,
			c1.g,
			c1.b,
		};
        temp.push_back(new_face1);

        return true;
    }
	
	return false;
}

enum LineType {
    Vertical,
    GradientLower,
    GradientUpper
};

LineType getLineType(double dx, double dy) {
    if (dx == 0) return Vertical;
    else if (abs(dy) < abs(dx)) return GradientLower;
    else return GradientUpper;
}

void render_vertical_line(Scene &scene, Vec3_helper first, Vec3_helper second, bool zero_to_one) {
	if(!zero_to_one) {
		Vec3_helper temp;
		
		temp.x = second.x;
		temp.y = second.y;
		temp.z = second.z;
		temp.r = second.r;
		temp.g = second.g;
		temp.b = second.b;

		second.x = first.x;
		second.y = first.y;
		second.z = first.z;
		second.r = first.r;
		second.g = first.g;
		second.b = first.b;

		first.x = temp.x;
		first.y = temp.y;
		first.z = temp.z;
		first.r = temp.r;
		first.g = temp.g;
		first.b = temp.b;
	}

	double delta_z = (second.z - first.z) / static_cast<double>(second.y - first.y);
	double delta_red = (second.r - first.r) / static_cast<double>(second.y - first.y);
	double delta_green = (second.g - first.g) / static_cast<double>(second.y - first.y);
	double delta_blue = (second.b - first.b) / static_cast<double>(second.y - first.y);

	int y = first.y;
	while (y <= second.y) {
		if (first.z < depth[first.x][y]) {
			scene.image[first.x][y].r = first.r;
			scene.image[first.x][y].g = first.g;
			scene.image[first.x][y].b = first.b;
			depth[first.x][y] = first.z;
		}
		first.r += delta_red; 
		first.g += delta_green; 
		first.b += delta_blue; 
		first.z += delta_z; 
		y++; 
	}
	
}

void render_lower_gradient(Scene &scene, Vec3_helper first, Vec3_helper second, bool zero_to_one) {
	if(!zero_to_one) {
		Vec3_helper temp;
		
		temp.x = second.x;
		temp.y = second.y;
		temp.z = second.z;
		temp.r = second.r;
		temp.g = second.g;
		temp.b = second.b;

		second.x = first.x;
		second.y = first.y;
		second.z = first.z;
		second.r = first.r;
		second.g = first.g;
		second.b = first.b;

		first.x = temp.x;
		first.y = temp.y;
		first.z = temp.z;
		first.r = temp.r;
		first.g = temp.g;
		first.b = temp.b;
	}	
	

	
	int delta_x = second.x - first.x;

	int y = first.y;
	int delta_y = first.y - second.y;
	int change_in_y = (delta_y > 0) ? -1 : 1;
	delta_y = (delta_y > 0) ? -delta_y : delta_y;

	int d = 2 * delta_y + delta_x;

	double z = first.z;
	double delta_z = (second.z - first.z) / static_cast<double>(delta_x);

	double delta_red = (second.r - first.r) / static_cast<double>(delta_x);
	double delta_green = (second.g - first.g) / static_cast<double>(delta_x);
	double delta_blue = (second.b - first.b) / static_cast<double>(delta_x);

	int x = first.x;
	while (x <= second.x) {
		if (z < depth[x][y]) {
			scene.image[x][y].r = first.r;
			scene.image[x][y].g = first.g;
			scene.image[x][y].b = first.b;
			depth[x][y] = z;
		}
		first.r += delta_red;
		first.g += delta_green;
		first.b += delta_blue; 

		d < 0 ? (y += change_in_y, d += 2 * (delta_y + delta_x)) : (d += 2 * delta_y);
		z += delta_z;
		x++; 
	}

}

void render_upper_gradient(Scene &scene, Vec3_helper first, Vec3_helper second, bool zero_to_one) {
	if(!zero_to_one) {
		Vec3_helper temp;
		
		temp.x = second.x;
		temp.y = second.y;
		temp.z = second.z;
		temp.r = second.r;
		temp.g = second.g;
		temp.b = second.b;

		second.x = first.x;
		second.y = first.y;
		second.z = first.z;
		second.r = first.r;
		second.g = first.g;
		second.b = first.b;

		first.x = temp.x;
		first.y = temp.y;
		first.z = temp.z;
		first.r = temp.r;
		first.g = temp.g;
		first.b = temp.b;
	}	
	
	int delta_y = second.y - first.y;

	int x = first.x;
	int delta_x = first.x - second.x;
	int change_in_x = (delta_x > 0) ? -1 : 1;
	delta_x = (delta_x > 0) ? -delta_x : delta_x;

	int d = 2 * delta_x + delta_y;

	double z = first.z;
	double delta_z = (second.z - first.z) / static_cast<double>(delta_y);

	double delta_red = (second.r - first.r) / static_cast<double>(delta_y);
	double delta_green = (second.g - first.g) / static_cast<double>(delta_y);
	double delta_blue = (second.b - first.b) / static_cast<double>(delta_y);

	int y = first.y;
	while (y <= second.y) {
		if (z < depth[x][y]) {
			scene.image[x][y].r = first.r;
			scene.image[x][y].g = first.g;
			scene.image[x][y].b = first.b;

			depth[x][y] = z;
		}

		first.r += delta_red;
		first.g += delta_green;
		first.b += delta_blue;


		d < 0 ? (x += change_in_x, d += 2 * (delta_x + delta_y)) : (d += 2 * delta_x);
		z += delta_z; 
		y++; 
	}

}

void rasterize_wire(Scene &scene, Camera &camera, vector<Vec4> &faces, Matrix4 &viewport_transformation) {
	vector<vector<Vec4_helper>> lines;
	vector<Vec4_helper> temp;
	
	if (to_clip(scene, camera, faces[0], faces[1], temp)) lines.push_back(temp);
	if (to_clip(scene, camera, faces[1], faces[2], temp)) lines.push_back(temp);
	if (to_clip(scene, camera, faces[2], faces[0], temp)) lines.push_back(temp);
	
	
	int i = 0;
	while (i < lines.size()) {
		vector<Vec3_helper> line;
		
		int j = 0;
		while (j < 2) { 
			Vec3_helper currline3;

			
			double values[4];
			double total;

			int k = 0;
			while (k < 4)
			{
				total = 0;
				int l = 0;

				while (l < 4)
				{
					double n;
					if(l == 0) n = lines[i][j].x;
					else if(l == 1) n = lines[i][j].y;
					else if (l == 2) n = lines[i][j].z;
					else n = lines[i][j].t;


					total += viewport_transformation.values[k][l] * n;
					l++;
				}

				values[k] = total;
				k++;
			}

			currline3.x = values[0];
			currline3.y = values[1];
			currline3.z = values[2];		
			currline3.r = lines[i][j].r;
			currline3.g = lines[i][j].g;
			currline3.b = lines[i][j].b;

			line.push_back(currline3);
			
			j++;
		}
		
		double delta_x = line[1].x - line[0].x;
		double delta_y = line[1].y - line[0].y;
		bool zero_to_one;

		switch (getLineType(delta_x, delta_y)) {
			case Vertical:
				zero_to_one = (delta_y < 0) ? false : true;
				render_vertical_line(scene, line[0], line[1], zero_to_one);
				break;

			case GradientLower:
				zero_to_one = (delta_x < 0) ? false : true;
				render_lower_gradient(scene, line[0], line[1], zero_to_one);
				break;

			case GradientUpper:
				zero_to_one = (delta_y < 0) ? false : true;
				render_upper_gradient(scene, line[0], line[1], zero_to_one);
				break;
		}

		i++;
	}

}

void Scene::forwardRenderingPipeline(Camera *camera)
{
	// TODO: Implement this function
	Matrix4 projection_transformation = getIdentityMatrix();
	camera_projection_transformation(*camera, projection_transformation);

	Matrix4 viewport_transformation = getIdentityMatrix();
	viewport_transformation.values[0][0] = static_cast<double>(camera->horRes) / 2;
	viewport_transformation.values[0][3] = static_cast<double>(camera->horRes - 1) / 2;
	viewport_transformation.values[1][1] = static_cast<double>(camera->verRes) / 2;
	viewport_transformation.values[1][3] = static_cast<double>(camera->verRes - 1) / 2;
	viewport_transformation.values[2][2] = 0.5;
	viewport_transformation.values[2][3] = 0.5;

	
	for(auto mesh: meshes) {

		Matrix4 model_transformation = getIdentityMatrix();
		find_model_transformation(*this, *mesh, model_transformation);
		model_transformation = multiplyMatrixWithMatrix(projection_transformation, model_transformation);

		
		for(auto triangle: mesh->triangles) {

			Vec3 vertices[3];
			vector<Vec4> faces;
			int i = 0;

			while(i < 3) {
				vertices[i] = *this->vertices[triangle.vertexIds[i] - 1];
				i++;
			}

			
			i = 0;
			while(i < 3) {
				Vec4 temp;
				temp.x = vertices[i].x;
				temp.y = vertices[i].y;
				temp.z = vertices[i].z;
				temp.t = 1;
				temp.colorId = vertices[i].colorId;

				temp = multiplyMatrixWithVec4(model_transformation, temp);

				if(camera->projectionType) {
					temp.x /= temp.t;
					temp.y /= temp.t;
					temp.z /= temp.t;
					temp.t = 1;
				}

				faces.push_back(temp);
				i++;
			}

			
			bool back_face = false;
			Vec3 v[3], temp1, temp2, temp3, temp4;
			i = 0;
			
			while(i < 3) {
				v[i].x = faces[i].x;
				v[i].y = faces[i].y;
				v[i].z = faces[i].z;
				v[i].colorId = faces[i].colorId; 
				i++;
			}

			temp1 = subtractVec3(v[1], v[0]);
			temp2 = subtractVec3(v[2], v[0]);
			temp1 = crossProductVec3(temp1, temp2);

			Vec3 dot;
			dot.x = dot.y = 0;
			dot.z = 1;


			if(dotProductVec3(dot, temp1) < 0) back_face = true;
			else back_face = false;

			if(this -> cullingEnabled && back_face) continue;

			if(mesh -> type) rasterize_mesh(*this, *camera, faces, viewport_transformation);
			else rasterize_wire(*this, *camera, faces, viewport_transformation);

			
		}

	}

}