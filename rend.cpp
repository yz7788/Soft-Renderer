/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

#define PI (float) 3.14159265358979323846

/* Vector Helper Function */
float GetNorm(const GzCoord& v)
{
	return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

float DotProduct(const GzCoord& u, const GzCoord& v)
{
	return u[0] * v[0] + u[1] * v[1] + u[2] * v[2];
}

void CrossProduct(const GzCoord& u, const GzCoord& v, GzCoord* result)
{
	(*result)[0] = u[1] * v[2] - u[2] * v[1];
	(*result)[1] = u[2] * v[0] - u[0] * v[2];
	(*result)[2] = u[0] * v[1] - u[1] * v[0];
	return;
}

void PlusVector(const GzCoord& v1, const GzCoord& v2, GzCoord* result)
{
	for (int i = 0; i < 3; i++)
	{
		(*result)[i] = v1[i] + v2[i];
	}
	return;
}

void MinusVector(const GzCoord& v1, const GzCoord& v2, GzCoord* result)
{
	for (int i = 0; i < 3; i++)
	{
		(*result)[i] = v1[i] - v2[i];
	}
	return;
}

void MultiplyNum(const GzCoord& v, const float& n, GzCoord* result)
{
	for (int i = 0; i < 3; i++)
	{
		(*result)[i] = v[i] * n;
	}
	return;
}

int DivideNum(const GzCoord& v, const float& n, GzCoord* result)
{
	if (n != 0)
	{
		for (int i = 0; i < 3; i++)
		{
			(*result)[i] = v[i] / n;
		}
		return GZ_SUCCESS;
	}
	else
		return GZ_FAILURE;
}

void Transform(const GzMatrix& m, GzCoord* v)
{
	float result[4];
	for (int i = 0; i < 4; i++)
	{
		result[i] = m[i][0] * (*v)[0] + m[i][1] * (*v)[1] + m[i][2] * (*v)[2] + m[i][3];
	}
	for (int i = 0; i < 3; i++)
	{
		(*v)[i] = result[i] / result[3];
	}
	return;
}

void Transform(const GzMatrix& m, const GzCoord& v, GzCoord* result)
{
	float temp[4];
	for (int i = 0; i < 4; i++)
	{
		temp[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2] + m[i][3];
	}
	for (int i = 0; i < 3; i++)
	{
		(*result)[i] = temp[i] / temp[3];
	}
	return;
}

void Normalize(GzCoord* v)
{
	float norm = GetNorm(*v);

	(*v)[0] = (*v)[0] / norm;
	(*v)[1] = (*v)[1] / norm;
	(*v)[2] = (*v)[2] / norm;

	return;
}

void MultiplyMatrix(const GzMatrix& m1, const GzMatrix& m2, GzMatrix* result)
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			(*result)[i][j] = 0;
			for (int k = 0; k < 4; k++)
			{
				(*result)[i][j] += m1[i][k] * m2[k][j];
			}
		}
	}
}


int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/* HW 3.1
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = 0;
		}
	}

	mat[0][0] = mat[3][3] = 1;
	mat[1][1] = mat[2][2] = cos(degree * (PI / 180));
	mat[1][2] = -sin(degree * (PI / 180));
	mat[2][1] = sin(degree * (PI / 180));

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/* HW 3.2
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = 0;
		}
	}

	mat[1][1] = mat[3][3] = 1;
	mat[0][0] = mat[2][2] = cos(degree * (PI / 180));
	mat[0][2] = sin(degree * (PI / 180));
	mat[2][0] = -sin(degree * (PI / 180));

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/* HW 3.3
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			mat[i][j] = 0;
		}
	}

	mat[2][2] = mat[3][3] = 1;
	mat[0][0] = mat[1][1] = cos(degree * (PI / 180));
	mat[1][0] = sin(degree * (PI / 180));
	mat[0][1] = -sin(degree * (PI / 180));

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/* HW 3.4
	// Create translation matrix
	// Pass back the matrix using mat value
	*/

	mat[0][0] = 1; mat[0][1] = 0; mat[0][2] = 0; mat[0][3] = translate[0];
	mat[1][0] = 0; mat[1][1] = 1; mat[1][2] = 0; mat[1][3] = translate[1];
	mat[2][0] = 0; mat[2][1] = 0; mat[2][2] = 1; mat[2][3] = translate[2];
	mat[3][0] = 0; mat[3][1] = 0; mat[3][2] = 0; mat[3][3] = 1;

	return GZ_SUCCESS;
}


int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/* HW 3.5
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/

	mat[0][0] = scale[0];  mat[0][1] = 0;         mat[0][2] = 0;         mat[0][3] = 0;
	mat[1][0] = 0;         mat[1][1] = scale[1];  mat[1][2] = 0;         mat[1][3] = 0;
	mat[2][0] = 0;         mat[2][1] = 0;         mat[2][2] = scale[2];  mat[2][3] = 0;
	mat[3][0] = 0;         mat[3][1] = 0;         mat[3][2] = 0;         mat[3][3] = 1;

	return GZ_SUCCESS;
}



GzRender::GzRender(int xRes, int yRes)
{
	/* HW1.1 create a framebuffer for MS Windows display:
 -- set display resolution
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- allocate memory for pixel buffer
 */

	this->xres = xRes;
	this->yres = yRes;
	this->framebuffer = (char*)malloc(3 * xRes * yRes * sizeof(char));
	this->pixelbuffer = (GzPixel*)malloc(xRes * yRes * sizeof(GzPixel));


	/* HW 3.6
	- setup Xsp and anything only done once
	- init default camera
	*/

	this->matlevel = -1;

	this->numlights = 0;

	this->m_camera.FOV = DEFAULT_FOV * (PI / 180);

	this->m_camera.position[0] = DEFAULT_IM_X;
	this->m_camera.position[1] = DEFAULT_IM_Y;
	this->m_camera.position[2] = DEFAULT_IM_Z;

	this->m_camera.lookat[0] = 0;
	this->m_camera.lookat[1] = 0;
	this->m_camera.lookat[2] = 0;

	this->m_camera.worldup[0] = 0;
	this->m_camera.worldup[1] = 1;
	this->m_camera.worldup[2] = 0;

	//Set up Xsp
	Xsp[0][0] = xres / 2;	Xsp[0][1] = 0.0f;		Xsp[0][2] = 0.0f;	 Xsp[0][3] = xres / 2;
	Xsp[1][0] = 0.0f;		Xsp[1][1] = -yres / 2;	Xsp[1][2] = 0.0f;	 Xsp[1][3] = yres / 2;
	Xsp[2][0] = 0.0f;		Xsp[2][1] = 0.0f;		Xsp[2][2] = MAXINT;	 Xsp[2][3] = 0.0f;
	Xsp[3][0] = 0.0f;		Xsp[3][1] = 0.0f;		Xsp[3][2] = 0.0f;	 Xsp[3][3] = 1.0f;
}

GzRender::~GzRender()
{
	/* HW1.2 clean up, free buffer memory */

	if (pixelbuffer != nullptr)
	{
		delete pixelbuffer;
	}
	if (framebuffer != nullptr)
	{
		delete framebuffer;
	}
}

int GzRender::GzDefault()
{
	/* HW1.3 set pixel buffer to some default values - start a new frame */

	for (int i = 0; i < xres; i++)
	{
		for (int j = 0; j < yres; j++)
		{
			int idx = ARRAY(i, j);
			GzPixel* currPixel = &pixelbuffer[idx];

			currPixel->red = 4000;
			currPixel->green = 3424;
			currPixel->blue = 3072;
			currPixel->alpha = 1;
			currPixel->z = INT_MAX;
		}
	}
	return GZ_SUCCESS;
}


int GzRender::GzBeginRender()
{
	/* HW 3.7
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	GzDefault();

	/* Set up Xpi */
	m_camera.Xpi[0][0] = 1.0f;  m_camera.Xpi[0][1] = 0.0f;  m_camera.Xpi[0][2] = 0.0f;									    m_camera.Xpi[0][3] = 0.0f;
	m_camera.Xpi[1][0] = 0.0f;  m_camera.Xpi[1][1] = 1.0f;  m_camera.Xpi[1][2] = 0.0f;									    m_camera.Xpi[1][3] = 0.0f;
	m_camera.Xpi[2][0] = 0.0f;  m_camera.Xpi[2][1] = 0.0f;  m_camera.Xpi[2][2] = tan((m_camera.FOV / 2.0f) * PI / 180.0f);  m_camera.Xpi[2][3] = 0.0f;
	m_camera.Xpi[3][0] = 0.0f;  m_camera.Xpi[3][1] = 0.0f;  m_camera.Xpi[3][2] = tan((m_camera.FOV / 2.0f) * PI / 180.0f);  m_camera.Xpi[3][3] = 1.0f;

	/* Set up Xiw */
	GzCoord xVec, yVec, zVec;

	// Compute zVec
	GzCoord cl;

	MinusVector(m_camera.lookat, m_camera.position, &cl);

	float clNorm = GetNorm(cl);

	DivideNum(cl, clNorm, &zVec);

	// Compute yVec
	GzCoord up;

	float dotProductUpZ = DotProduct(m_camera.worldup, zVec);

	up[0] = m_camera.worldup[0] - (dotProductUpZ * zVec[0]);
	up[1] = m_camera.worldup[1] - (dotProductUpZ * zVec[1]);
	up[2] = m_camera.worldup[2] - (dotProductUpZ * zVec[2]);

	float modUp = GetNorm(up);

	DivideNum(up, modUp, &yVec);

	// Compute xVec
	CrossProduct(yVec, zVec, &xVec);

	m_camera.Xiw[0][0] = xVec[0];  m_camera.Xiw[0][1] = xVec[1];  m_camera.Xiw[0][2] = xVec[2];  m_camera.Xiw[0][3] = -DotProduct(xVec, m_camera.position);
	m_camera.Xiw[1][0] = yVec[0];  m_camera.Xiw[1][1] = yVec[1];  m_camera.Xiw[1][2] = yVec[2];  m_camera.Xiw[1][3] = -DotProduct(yVec, m_camera.position);
	m_camera.Xiw[2][0] = zVec[0];  m_camera.Xiw[2][1] = zVec[1];  m_camera.Xiw[2][2] = zVec[2];  m_camera.Xiw[2][3] = -DotProduct(zVec, m_camera.position);
	m_camera.Xiw[3][0] = 0;        m_camera.Xiw[3][1] = 0;        m_camera.Xiw[3][2] = 0;        m_camera.Xiw[3][3] = 1;

	/* Set up Xnorm (Xim) & Ximage (Xsm) */
	GzPushMatrix(Xsp, false);
	GzPushMatrix(m_camera.Xpi, false);
	GzPushMatrix(m_camera.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/* HW 3.8
	/*- overwrite renderer camera structure with new camera definition
	*/

	m_camera = camera;
	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(const GzMatrix& matrix, bool pushToXnorm)
{
	/* HW 3.9
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/

	if (matlevel >= MATLEVELS) return GZ_FAILURE;

	GzMatrix identity{ {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };

	matlevel++;
	if (matlevel == 0)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Ximage[matlevel][i][j] = matrix[i][j];
			}
		}
	}
	else
	{
		MultiplyMatrix(Ximage[matlevel - 1], matrix, &Ximage[matlevel]);
	}

	if (pushToXnorm)
	{
		GzMatrix unitaryMatrix{ {matrix[0][0], matrix[0][1], matrix[0][2], 0},
			{matrix[1][0], matrix[1][1], matrix[1][2], 0},
			{matrix[2][0], matrix[2][1], matrix[2][2], 0},
			{matrix[3][0], matrix[3][1], matrix[3][2], 1} };

		//normalize unitaryMatrix
		float norm = 1 / sqrt(matrix[0][0] * matrix[0][0] + matrix[0][1] * matrix[0][1] + matrix[0][2] * matrix[0][2]);

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				unitaryMatrix[i][j] *= norm;
			}
		}

		if (matlevel == 0)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Xnorm[matlevel][i][j] = unitaryMatrix[i][j];
				}
			}
		}
		else
		{
			MultiplyMatrix(Xnorm[matlevel - 1], unitaryMatrix, &Xnorm[matlevel]);
		}
	}
	else
	{
		if (matlevel == 0)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Xnorm[matlevel][i][j] = identity[i][j];
				}
			}
		}
		else
		{
			MultiplyMatrix(Xnorm[matlevel - 1], identity, &Xnorm[matlevel]);
		}
	}

	return GZ_SUCCESS;
}

int GzRender::GzPopMatrix()
{
	/* HW 3.10
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/

	if (matlevel < 0) return GZ_FAILURE;
	else {
		matlevel--;
	}

	return GZ_SUCCESS;
}

/*
Description:
This function is used to write pixel values into the pixel buffer;
Input:
@ int i, int j: x, y position;
@ GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a: red, green, blue, alpha value;
@ GzDepth z: depth value;
Output:
@ int returnValue: status;
*/
int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.4 write pixel values into the buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) 
	{
		int idx = ARRAY(i, j);
		GzPixel* currPixel = &pixelbuffer[idx];
		if (z < currPixel->z)
		{
			currPixel->red = r < 4095 ? r : 4095;
			currPixel->green = g < 4095 ? g : 4095;
			currPixel->blue = b < 4095 ? b : 4095;
			currPixel->alpha = a;
			currPixel->z = z;
		}
	}

	return GZ_SUCCESS;
}


int GzRender::GzGet(int i, int j, GzIntensity* r, GzIntensity* g, GzIntensity* b, GzIntensity* a, GzDepth* z)
{
	/* HW1.5 retrieve a pixel information from the pixel buffer */

	if (i >= 0 && i < xres && j >= 0 && j < yres) {
		int idx = ARRAY(i, j);
		GzPixel* currPixel = &pixelbuffer[idx];

		r = &currPixel->red;
		g = &currPixel->green;
		b = &currPixel->blue;
		a = &currPixel->alpha;
		z = &currPixel->z;
	}
	return GZ_SUCCESS;
}


int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* HW1.6 write image to ppm file -- "P6 %d %d 255\r" */

	fprintf(outfile, "P6 %d %d 255\r", xres, yres);
	for (int j = 0; j < yres; j++) {
		for (int i = 0; i < xres; i++) {
			int idx = ARRAY(i, j);
			GzPixel* currPixel = &pixelbuffer[idx];

			fprintf(outfile, "%c%c%c", currPixel->red >> 4, currPixel->green >> 4, currPixel->blue >> 4);
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* HW1.7 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int k = 0;
	for (int j = 0; j < yres; j++) {
		for (int i = 0; i < xres; i++) {
			int idx = ARRAY(i, j);
			GzPixel* currPixel = &pixelbuffer[idx];

			framebuffer[k++] = currPixel->blue >> 4;
			framebuffer[k++] = currPixel->green >> 4;
			framebuffer[k++] = currPixel->red >> 4;
		}
	}
	return GZ_SUCCESS;
}

/***********************************************/
/* HW2 methods: implement from here */

/*
Description:
This function is used to set renderer attribute states;
Input:
@ void numberAttributes: the number of attributes;
@ GzToken *nameList: a list of token, where the token indicates the rendering classification;
@ GzPointer *valueList: a list of attribute value;
Output:
@ int returnValue: status;
*/
int GzRender::GzPutAttribute(int numAttributes, GzToken* nameList, GzPointer* valueList)
{
	/* HW 2.1
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	for (int i = 0; i < numAttributes; i++)
	{
		switch (nameList[i])
		{
		case GZ_RGB_COLOR: 
		{
			GzColor* color = (GzColor*)(valueList[i]);

			this->flatcolor[0] = ctoi(*color[0]);
			this->flatcolor[1] = ctoi(*color[1]);
			this->flatcolor[2] = ctoi(*color[2]);

		} break;
		case GZ_DIRECTIONAL_LIGHT:
		{
			this->lights[numlights].direction[0] = ((GzLight*)(valueList[i]))->direction[0];
			this->lights[numlights].direction[1] = ((GzLight*)(valueList[i]))->direction[1];
			this->lights[numlights].direction[2] = ((GzLight*)(valueList[i]))->direction[2];

			this->lights[numlights].color[0] = ((GzLight*)(valueList[i]))->color[0];
			this->lights[numlights].color[1] = ((GzLight*)(valueList[i]))->color[1];
			this->lights[numlights].color[2] = ((GzLight*)(valueList[i]))->color[2];

			numlights++;
		} break;
		case GZ_AMBIENT_LIGHT:
		{
			this->ambientlight.direction[0] = ((GzLight*)(valueList[i]))->direction[0];
			this->ambientlight.direction[1] = ((GzLight*)(valueList[i]))->direction[1];
			this->ambientlight.direction[2] = ((GzLight*)(valueList[i]))->direction[2];

			this->ambientlight.color[0] = ((GzLight*)(valueList[i]))->color[0];
			this->ambientlight.color[1] = ((GzLight*)(valueList[i]))->color[1];
			this->ambientlight.color[2] = ((GzLight*)(valueList[i]))->color[2];

		} break;
		case GZ_INTERPOLATE:
		{
			this->interp_mode = *(int*)(valueList[i]);

		} break;
		case GZ_AMBIENT_COEFFICIENT:
		{
			this->Ka[0] = (*(GzColor*)(valueList[i]))[0];
			this->Ka[1] = (*(GzColor*)(valueList[i]))[1];
			this->Ka[2] = (*(GzColor*)(valueList[i]))[2];

		} break;
		case GZ_DIFFUSE_COEFFICIENT:
		{
			this->Kd[0] = (*(GzColor*)(valueList[i]))[0];
			this->Kd[1] = (*(GzColor*)(valueList[i]))[1];
			this->Kd[2] = (*(GzColor*)(valueList[i]))[2];

		} break;
		case GZ_SPECULAR_COEFFICIENT:
		{
			this->Ks[0] = (*(GzColor*)(valueList[i]))[0];
			this->Ks[1] = (*(GzColor*)(valueList[i]))[1];
			this->Ks[2] = (*(GzColor*)(valueList[i]))[2];

		} break;
		case GZ_DISTRIBUTION_COEFFICIENT:
		{
			this->spec = *(float*)(valueList[i]);
		} break;
		}
	}
	return GZ_SUCCESS;

}

#pragma region PutTriangle

float GetA(const GzCoord& v1, const GzCoord& v2)
{
	return v2[1] - v1[1];
}

float GetB(const GzCoord& v1, const GzCoord& v2)
{
	return -(v2[0] - v1[0]);
}

float GetC(const GzCoord& v1, const GzCoord& v2)
{
	return ((v2[0] - v1[0]) * v1[1] - (v2[1] - v1[1]) * v1[0]);
}

bool PointInsideLine(float x1, float y1, float x2, float y2, float x, float y)
{
	float A = y2 - y1;
	float B = (-(x2 - x1));
	float C = ((x2 - x1) * y1 - (y2 - y1) * x1);

	return ((A * x + B * y + C) < 0);
}

/*
Description:
This function is used to get an interpolated Z value of a specific point, given the corresponding x,y value and the vertices information of the triangle
Input:
@ GzCoord** v1, GzCoord** v2, GzCoord** v3: 3 triangle vertices
@ int i, int j: x, y value of the point
Output:
@ float returnValue: interpolated z-value
*/
float InterpolateZ(GzCoord* v1, GzCoord* v2, GzCoord* v3, const int& i, const int& j)
{
	float A = ((*v2)[1] - (*v1)[1]) * ((*v3)[2] - (*v1)[2]) - ((*v2)[2] - (*v1)[2]) * ((*v3)[1] - (*v1)[1]);
	float B = ((*v2)[2] - (*v1)[2]) * ((*v3)[0] - (*v1)[0]) - ((*v2)[0] - (*v1)[0]) * ((*v3)[2] - (*v1)[2]);
	float C = ((*v2)[0] - (*v1)[0]) * ((*v3)[1] - (*v1)[1]) - ((*v2)[1] - (*v1)[1]) * ((*v3)[0] - (*v1)[0]);
	float D = -(A * (*v1)[0] + B * (*v1)[1] + C * (*v1)[2]);

	return (-D - A * (float)i - B * (float)j) / C;
}

void InterpolateColor(const GzCoord& v1, const GzCoord& v2, const GzCoord& v3, const GzColor& c1, const GzColor& c2, const GzColor& c3, const int& i, const int& j, GzColor* c)
{
	// calculate R
	GzCoord R1{ v1[0], v1[1], c1[0] };
	GzCoord R2{ v2[0], v2[1], c2[0] };
	GzCoord R3{ v3[0], v3[1], c3[0] };

	(*c)[0] = InterpolateZ(&R1, &R2, &R3, i, j);

	// calculate G
	GzCoord G1{ v1[0], v1[1], c1[1] };
	GzCoord G2{ v2[0], v2[1], c2[1] };
	GzCoord G3{ v3[0], v3[1], c3[1] };

	(*c)[1] = InterpolateZ(&G1, &G2, &G3, i, j);

	// calculate B
	GzCoord B1{ v1[0], v1[1], c1[2] };
	GzCoord B2{ v2[0], v2[1], c2[2] };
	GzCoord B3{ v3[0], v3[1], c3[2] };

	(*c)[2] = InterpolateZ(&B1, &B2, &B3, i, j);

	return;
}


/*
Description:
This function is used to calculate color based on color equation including specular color, diffuse color, and ambient color;
Input:
@ int numlights, GzLight* lights, GzLight* ambientlight, GzColor Ka, GzColor Kd, GzColor Ks, float spec, GzCoord n;
Output:
@ GzColor c1, GzColor c2, GzColor c3;
*/
void GzRender::CalculateColor(const GzCoord& n, GzColor* c)
{
	//ambient component
	
	GzCoord E{ 0, 0, -1 };

	GzColor ambientComponent{ this->Ka[0] * this->ambientlight.color[0] , this->Ka[1] * this->ambientlight.color[1] , this->Ka[2] * this->ambientlight.color[2] };
	GzColor diffuseComponent{0, 0, 0};
	GzColor specularComponent{0, 0, 0};

	float nDotE = DotProduct(n, E);
	for (int i = 0; i < numlights; i++)
	{
		// diffuse component 
		GzCoord L{ this->lights[i].direction[0], this->lights[i].direction[1], this->lights[i].direction[2]};
		Normalize(&L);

		float nDotL = DotProduct(n, L);
		
		//specular component 
		GzCoord R{ 0, 0, 0 };

		if (nDotL > 0 && nDotE > 0)
		{
			//R = normalized(2 * n * dotproduct(L, n) - L);
			MultiplyNum(n, 2 * DotProduct(L, n), &R);
			MinusVector(R, L, &R);
			Normalize(&R);
		}
		else if (nDotL < 0 && nDotE < 0)
		{
			GzCoord flipN{ -n[0], -n[1], -n[2] };

			nDotL = DotProduct(n, L);
			
			MultiplyNum(n, 2 * DotProduct(L, flipN), &R);
			MinusVector(R, L, &R);
			Normalize(&R);
		}
		
		float d = max(0.0f, DotProduct(n, L));

		diffuseComponent[0] += this->Kd[0] * this->lights[i].color[0] * d;
		diffuseComponent[1] += this->Kd[1] * this->lights[i].color[1] * d;
		diffuseComponent[2] += this->Kd[2] * this->lights[i].color[2] * d;

		float s = min(1.0f, max(0.0f, DotProduct(R, E)));

		specularComponent[0] += this->Ks[0] * this->lights[i].color[0] * pow(s, this->spec);
		specularComponent[1] += this->Ks[1] * this->lights[i].color[1] * pow(s, this->spec);
		specularComponent[2] += this->Ks[2] * this->lights[i].color[2] * pow(s, this->spec);
	}
	
	(*c)[0] = ambientComponent[0] + diffuseComponent[0] + specularComponent[0];
	(*c)[1] = ambientComponent[1] + diffuseComponent[1] + specularComponent[1];
	(*c)[2] = ambientComponent[2] + diffuseComponent[2] + specularComponent[2];

	return;
}

int GzRender::GzPutTriangle(int	numParts, GzToken* nameList, GzPointer* valueList)
/* numParts - how many names and values */
{
	/*
	-- Pass in a triangle description with tokens and values corresponding to
		  GZ_NULL_TOKEN:		do nothing - no values
		  GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/

	GzCoord* v1;
	GzCoord* v2;
	GzCoord* v3;

	GzCoord* n1;
	GzCoord* n2;
	GzCoord* n3;

	for (int i = 0; i < numParts; i++)
	{
		switch (nameList[i])
		{
		case GZ_NULL_TOKEN:
			break;
		case GZ_POSITION: {

			GzCoord* vertices = (GzCoord*)valueList[i];
			v1 = vertices;
			v2 = vertices + 1;
			v3 = vertices + 2;

		} break;
		case GZ_NORMAL: {

			GzCoord* normals = (GzCoord*)valueList[i];
			n1 = normals;
			n2 = normals + 1;
			n3 = normals + 2;
		
			//std::string s = std::to_string((*n1)[2]);
			//CString ss(s.c_str());
			//AfxMessageBox(ss);
		} 
		}
	}

	//transform v1,v2,v3 into screen space 
	Transform(Ximage[matlevel], v1);
	Transform(Ximage[matlevel], v2);
	Transform(Ximage[matlevel], v3);

	VertexSorting(&v1, &v2, &v3, &n1, &n2, &n3);

	// transform n1, n2, n3 into image space
	Transform(Xnorm[matlevel], n1);
	Transform(Xnorm[matlevel], n2);
	Transform(Xnorm[matlevel], n3);

	Normalize(n1);
	Normalize(n2);
	Normalize(n3);

	// calculate color for 3 vertices
	GzColor c1{ 0, 0, 0 };
	GzColor c2{ 0, 0, 0 };
	GzColor c3{ 0, 0, 0 };

	CalculateColor(*n1, &c1);
	CalculateColor(*n2, &c2);
	CalculateColor(*n3, &c3);

	if (v1 && v2 && v3 && n1 && n2 && n3)
	{
		int x_min = min(min((*v1)[0], (*v2)[0]), (*v3)[0]);
		int y_min = (int)(*v1)[1];
		int x_max = max(max((*v1)[0], (*v2)[0]), (*v3)[0]);
		int y_max = max((*v2)[1], (*v3)[1]);

		for (int i = x_min; i <= x_max; i++)
		{
			for (int j = y_min; j <= y_max; j++)
			{
				if (PointInsideLine((*v1)[0], (*v1)[1], (*v2)[0], (*v2)[1], i, j)
					&& PointInsideLine((*v2)[0], (*v2)[1], (*v3)[0], (*v3)[1], i, j)
					&& PointInsideLine((*v3)[0], (*v3)[1], (*v1)[0], (*v1)[1], i, j))
				{
					if (this->interp_mode == GZ_FLAT) {
						//std::string s = std::to_string(InterpolateZ(v1, v2, v3, i, j));
						//CString ss(s.c_str());
						//AfxMessageBox(ss);
						this->flatcolor[0] = c1[0];
						this->flatcolor[1] = c1[1];
						this->flatcolor[2] = c1[2];
					}
					else if (this->interp_mode == GZ_COLOR)
					{
						InterpolateColor(*v1, *v2, *v3, c1, c2, c3, i, j, &(this->flatcolor));
					}
					else if (this->interp_mode == GZ_NORMALS)
					{
						GzCoord tempN{0, 0, 0};
						GzColor tempC{ 0, 0, 0 };
						InterpolateColor(*v1, *v2, *v3, *n1, *n2, *n3, i, j, &tempN);
						CalculateColor(tempN, &tempC);
						this->flatcolor[0] = tempC[0];
						this->flatcolor[1] = tempC[1];
						this->flatcolor[2] = tempC[2];
					}
					else
					{
						CString s = "Wrong Interpolation Mode";
						AfxMessageBox(s);
					}
				
					GzPut(i, j, ctoi(flatcolor[0]), ctoi(flatcolor[1]), ctoi(flatcolor[2]), 1, InterpolateZ(v1, v2, v3, i, j));
				}
			}
		}
	}
	else
	{
		CString s = "Not enough vertices information";
		AfxMessageBox(s);
		return GZ_FAILURE;
	}

	return GZ_SUCCESS;
}


/*

bool sortHelper(const GzCoord* vi, const GzCoord* vj)
{
	return ((*vi)[1] < (*vj)[1]);
}

void GzRender::VertexSorting(GzCoord** v1, GzCoord** v2, GzCoord** v3)
{
	// start by sorting vert Y-coords 
	std::vector<GzCoord*> v{ *v1, *v2, *v3 };
	std::sort(v.begin(), v.end(), sortHelper);

	*v1 = v[0];
	*v2 = v[1];
	*v3 = v[2];

	//find L/R relationship, sort there vertices 
	if ((*v[0])[1] < (*v[1])[1])
	{
		// special case where there are x-axis's parallel line
		if ((*v[0])[0] == (*v[1])[0])
		{
			if ((*v[2])[0] > (*v[0])[0])
			{
				*v1 = v[0];
				*v2 = v[2];
				*v3 = v[1];
			}
		}
		else if ((*v[0])[0] == (*v[2])[0])
		{
			if ((*v[1])[0] < (*v[0])[0])
			{
				*v1 = v[0];
				*v2 = v[2];
				*v3 = v[1];
			}
		}
		// use mid-Y vert v[1], find X-coord of point P which is the intersection of edge v[0]-v[2] & x-axis's parallel line across v[1]  
		// compare point P and v[1], the one has greater x-coordinate is the right edge
		else
		{

			float A = GetA(*v[2], *v[0]);
			float B = GetB(*v[2], *v[0]);
			float C = GetC(*v[2], *v[0]);
			float xp = (-C - B * ((*v[1])[1])) / A;

			if (xp > (*v[1])[0])
			{
				*v1 = v[0];
				*v2 = v[2];
				*v3 = v[1];
			}
		}
	}
	// special case where there are y-axis's parallel line
	else
	{
		if ((*v[0])[0] > (*v[1])[0])
		{
			*v1 = v[0];
			*v2 = v[2];
			*v3 = v[1];
		}
	}

	return;
}

*/

bool sortHelper(const std::pair<GzCoord*, GzCoord*> pi, const std::pair<GzCoord*, GzCoord*> pj)
{
	return ((*pi.first)[1] < (*pj.first)[1]);
}

/*
Description:
This function is used to sort 3 vertices of a triangle in closewise order
Input:
@ GzCoord** v1, GzCoord** v2, GzCoord** v3: 3 vertices of a triangle
Output:
@ GzCoord** v1, GzCoord** v2, GzCoord** v3: 3 vertices of a triangle
*/

void GzRender::VertexSorting(GzCoord** v1, GzCoord** v2, GzCoord** v3, GzCoord** n1, GzCoord** n2, GzCoord** n3)
{
	//start by sorting vert Y-coords 
	std::pair<GzCoord*, GzCoord*> p1{ *v1, *n1};
	std::pair<GzCoord*, GzCoord*> p2{ *v2, *n2};
	std::pair<GzCoord*, GzCoord*> p3{ *v3, *n3};

	std::vector<std::pair<GzCoord*, GzCoord*>> p{ p1, p2, p3 };
	std::sort(p.begin(), p.end(), sortHelper);

	p1 = p[0];
	p2 = p[1];
	p3 = p[2];

	//find L/R relationship, sort three vertices
	if ((*(p[0].first))[1] < (*(p[1].first))[1])
	{
		// special case where there are x-axis's parallel line
		if ((*(p[0].first))[0] == (*(p[1].first))[0])
		{
			if ((*(p[2].first))[0] > (*(p[0].first))[0])
			{
				p1 = p[0];
				p2 = p[2];
				p3 = p[1];
			}
		}

		else if ((*(p[0].first))[0] == (*(p[2].first))[0])
		{
			if ((*(p[1].first))[0] < (*(p[0].first))[0])
			{
				p1 = p[0];
				p2 = p[2];
				p3 = p[1];
			}
		}

		// use mid-Y vert v[1], find X-coord of point P which is the intersection of edge v[0]-v[2] & x-axis's parallel line across v[1]  
		// compare point P and v[1], the one has greater x-coordinate is the right edge
		else
		{
			float A = (*(p[0].first))[1] - (*(p[2].first))[1];
			float B = (*(p[2].first))[0] - (*(p[0].first))[0];
			float C = ((*(p[0].first))[0] - (*(p[2].first))[0]) * (*(p[2].first))[1] - ((*(p[0].first))[1] - (*(p[2].first))[1]) * (*(p[2].first))[0];

			float xp = (-C - B * (*(p[1].first))[1]) / A;

			if (xp > (*(p[1].first))[0])
			{
				p1 = p[0];
				p2 = p[2];
				p3 = p[1];
			}
		}
	}
	// special case where there are y-axis's parallel line
	else
	{
		if ((*(p[0].first))[0] > (*(p[1].first))[0])
		{
			p1 = p[0];
			p2 = p[2];
			p3 = p[1];
		}
	}

	*v1 = p1.first;
	*v2 = p2.first;
	*v3 = p3.first;
	*n1 = p1.second;
	*n2 = p2.second;
	*n3 = p3.second;

	return;
}

#pragma endregion