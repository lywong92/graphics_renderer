#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

#define PI (float) 3.14159265358979323846

static bool cull = false;

// debug function to print out matrix values
void PrintMatrix(GzMatrix matrix)
{
	TRACE("PrintMatrix:\n");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			TRACE("matrix[%d][%d]: %f\n", i, j, matrix[i][j]);
		}
	}
}

float dot(GzCoord v1, GzCoord v2)
{
	return v1[X] * v2[X] + v1[Y] * v2[Y] + v1[Z] * v2[Z];
}

int GzRender::GzRotXMat(float degree, GzMatrix mat)
{
	/*
	// Create rotate matrix : rotate along x axis
	// Pass back the matrix using mat value
	*/
	mat[1][1] = (float)cos((degree * PI) / 180);
	mat[1][2] = (float)(-sin((degree * PI) / 180));
	mat[2][1] = (float)sin((degree * PI) / 180);
	mat[2][2] = (float)cos((degree * PI) / 180);

	return GZ_SUCCESS;
}

int GzRender::GzRotYMat(float degree, GzMatrix mat)
{
	/*
	// Create rotate matrix : rotate along y axis
	// Pass back the matrix using mat value
	*/
	mat[0][0] = (float)cos((degree * PI) / 180);
	mat[0][2] = (float)sin((degree * PI) / 180);
	mat[2][0] = (float)(-sin((degree * PI) / 180));
	mat[2][2] = (float)cos((degree * PI) / 180);

	return GZ_SUCCESS;
}

int GzRender::GzRotZMat(float degree, GzMatrix mat)
{
	/*
	// Create rotate matrix : rotate along z axis
	// Pass back the matrix using mat value
	*/
	mat[0][0] = (float)cos((degree * PI) / 180);
	mat[0][1] = (float)(-sin((degree * PI) / 180));
	mat[1][0] = (float)sin((degree * PI) / 180);
	mat[1][1] = (float)cos((degree * PI) / 180);

	return GZ_SUCCESS;
}

int GzRender::GzTrxMat(GzCoord translate, GzMatrix mat)
{
	/*
	// Create translation matrix
	// Pass back the matrix using mat value
	*/
	mat[0][3] = translate[X];
	mat[1][3] = translate[Y];
	mat[2][3] = translate[Z];

	return GZ_SUCCESS;
}

int GzRender::GzScaleMat(GzCoord scale, GzMatrix mat)
{
	/*
	// Create scaling matrix
	// Pass back the matrix using mat value
	*/
	mat[0][0] = scale[X];
	mat[1][1] = scale[Y];
	mat[2][2] = scale[Z];

	return GZ_SUCCESS;
}


GzRender::GzRender(int xRes, int yRes)
{
	/* create a framebuffer for MS Windows display:
	-- set display resolution
	-- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
	-- allocate memory for pixel buffer
	*/

	xres = xRes;
	yres = yRes;

	pixelbuffer = new GzPixel[xRes*yRes];
	// framebuffer stores BGR only for ouput to ppm file
	framebuffer = new char[xRes*yRes * 3];

	/*
	- setup Xsp and anything only done once
	- init default camera
	*/
	// create Xsp based on resolution
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			Xsp[i][j] = 0.0;
		}
	}
	Xsp[0][0] = (float)(xRes / 2);
	Xsp[0][3] = (float)(xRes / 2);
	Xsp[1][1] = (float)(-yRes / 2);
	Xsp[1][3] = (float)(yRes / 2);
	Xsp[2][2] = (float)MAXINT;
	Xsp[3][3] = (float)1;

	// default camera values
	m_camera.position[X] = DEFAULT_IM_X;
	m_camera.position[Y] = DEFAULT_IM_Y;
	m_camera.position[Z] = DEFAULT_IM_Z;

	m_camera.lookat[X] = 0.0;
	m_camera.lookat[Y] = 0.0;
	m_camera.lookat[Z] = 0.0;

	m_camera.worldup[X] = 0.0;
	m_camera.worldup[Y] = 1.0;
	m_camera.worldup[Z] = 0.0;

	m_camera.FOV = DEFAULT_FOV;

	matlevel = 0;
	numlights = 0;
}

GzRender::~GzRender()
{
	/* clean up, free buffer memory */
	delete[] pixelbuffer;
	delete[] framebuffer;
}

int GzRender::GzDefault()
{
	/* set pixel buffer to some default values - start a new frame */

	// loop through every pixel and put data in pixelbuffer
	for (int i = 0; i < xres * yres; i++)
	{
		pixelbuffer[i].red = 1234;
		pixelbuffer[i].green = 1234;
		pixelbuffer[i].blue = 1234;
		pixelbuffer[i].alpha = 1;
		pixelbuffer[i].z = MAXINT;
	}

	return GZ_SUCCESS;
}

int GzRender::GzBeginRender()
{
	/*
	- setup for start of each frame - init frame buffer color,alpha,z
	- compute Xiw and projection xform Xpi from camera definition
	- init Ximage - put Xsp at base of stack, push on Xpi and Xiw
	- now stack contains Xsw and app can push model Xforms when needed
	*/
	GzMatrix identity =
	{
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, 1.0, 0.0,
		0.0, 0.0, 0.0, 1.0
	};

	// loop through every pixel and put data in pixelbuffer
	for (int i = 0; i < (xres*yres); i++)
	{
		pixelbuffer[i].red = 1234;
		pixelbuffer[i].green = 1234;
		pixelbuffer[i].blue = 1234;
		pixelbuffer[i].alpha = 1;
		pixelbuffer[i].z = MAXINT;
	}

	// Xpi
	float one_over_d = (float)(tan((m_camera.FOV * PI) / 360.0));
	GzMatrix Xpi_copy =
	{
		1.0, 0.0, 0.0, 0.0,
		0.0, 1.0, 0.0, 0.0,
		0.0, 0.0, one_over_d, 0.0,
		0.0, 0.0, one_over_d, 1.0
	};
	memcpy(m_camera.Xpi, Xpi_copy, sizeof(Xpi_copy));

	// Xiw
	GzCoord camera_z;
	GzCoord camera_y;
	GzCoord camera_x;
	float cl_length, up_prime_length, up_dot_z;

	// calculate z
	camera_z[X] = m_camera.lookat[X] - m_camera.position[X];
	camera_z[Y] = m_camera.lookat[Y] - m_camera.position[Y];
	camera_z[Z] = m_camera.lookat[Z] - m_camera.position[Z];
	cl_length = (float)(sqrt(pow(camera_z[X], 2.0) +
		pow(camera_z[Y], 2.0) + pow(camera_z[Z], 2.0)));
	camera_z[X] = camera_z[X] / cl_length;
	camera_z[Y] = camera_z[Y] / cl_length;
	camera_z[Z] = camera_z[Z] / cl_length;

	//calculate y
	up_dot_z = dot(m_camera.worldup, camera_z);
	camera_y[X] = m_camera.worldup[X] - (up_dot_z * camera_z[X]);
	camera_y[Y] = m_camera.worldup[Y] - (up_dot_z * camera_z[Y]);
	camera_y[Z] = m_camera.worldup[Z] - (up_dot_z * camera_z[Z]);
	up_prime_length = (float)(sqrt(pow(camera_y[X], 2.0) +
		pow(camera_y[Y], 2.0) + pow(camera_y[Z], 2.0)));
	camera_y[X] = camera_y[X] / up_prime_length;
	camera_y[Y] = camera_y[Y] / up_prime_length;
	camera_y[Z] = camera_y[Z] / up_prime_length;

	// calulate x
	camera_x[X] = (camera_y[Y] * camera_z[Z]) - (camera_y[Z] * camera_z[Y]);
	camera_x[Y] = (camera_y[Z] * camera_z[X]) - (camera_y[X] * camera_z[Z]);
	camera_x[Z] = (camera_y[X] * camera_z[Y]) - (camera_y[Y] * camera_z[X]);

	// create Xiw
	m_camera.Xiw[0][0] = camera_x[X];
	m_camera.Xiw[0][1] = camera_x[Y];
	m_camera.Xiw[0][2] = camera_x[Z];
	m_camera.Xiw[0][3] =
		((-camera_x[X]) * m_camera.position[X]) +
		((-camera_x[Y]) * m_camera.position[Y]) +
		((-camera_x[Z]) * m_camera.position[Z]);
	m_camera.Xiw[1][0] = camera_y[X];
	m_camera.Xiw[1][1] = camera_y[Y];
	m_camera.Xiw[1][2] = camera_y[Z];
	m_camera.Xiw[1][3] =
		((-camera_y[X]) * m_camera.position[X]) +
		((-camera_y[Y]) * m_camera.position[Y]) +
		((-camera_y[Z]) * m_camera.position[Z]);
	m_camera.Xiw[2][0] = camera_z[X];
	m_camera.Xiw[2][1] = camera_z[Y];
	m_camera.Xiw[2][2] = camera_z[Z];
	m_camera.Xiw[2][3] =
		((-camera_z[X]) * m_camera.position[X]) +
		((-camera_z[Y]) * m_camera.position[Y]) +
		((-camera_z[Z]) * m_camera.position[Z]);
	m_camera.Xiw[3][0] = 0;
	m_camera.Xiw[3][1] = 0;
	m_camera.Xiw[3][2] = 0;
	m_camera.Xiw[3][3] = 1;

	// push identity matrix to bottom of stack first
	memcpy(Ximage[0], identity, sizeof(identity));
	memcpy(Xnorm[0], identity, sizeof(identity));

	// push Xsp, Xpi, Xiw
	GzPushMatrix(Xsp);
	GzPushMatrix(m_camera.Xpi);
	GzPushMatrix(m_camera.Xiw);

	return GZ_SUCCESS;
}

int GzRender::GzPutCamera(GzCamera camera)
{
	/*
	/*- overwrite renderer camera structure with new camera definition
	*/
	m_camera.position[X] = camera.position[X];
	m_camera.position[Y] = camera.position[Y];
	m_camera.position[Z] = camera.position[Z];

	m_camera.lookat[X] = camera.lookat[X];
	m_camera.lookat[Y] = camera.lookat[Y];
	m_camera.lookat[Z] = camera.lookat[Z];

	m_camera.worldup[X] = camera.worldup[X];
	m_camera.worldup[Y] = camera.worldup[Y];
	m_camera.worldup[Z] = camera.worldup[Z];

	m_camera.FOV = camera.FOV;

	return GZ_SUCCESS;
}

int GzRender::GzPushMatrix(GzMatrix	matrix)
{
	/*
	- push a matrix onto the Ximage stack
	- check for stack overflow
	*/
	float scale;

	if (matlevel < MATLEVELS)
	{
		// Ximage
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				Ximage[matlevel + 1][i][j] =
					Ximage[matlevel][i][0] * matrix[0][j]
					+ Ximage[matlevel][i][1] * matrix[1][j]
					+ Ximage[matlevel][i][2] * matrix[2][j]
					+ Ximage[matlevel][i][3] * matrix[3][j];
			}
		}

		// Xnorm
		// if Xsp or Xpi, dont push
		if (matrix[2][2] == MAXINT || matrix[3][2] != 0.0)
		{
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Xnorm[matlevel + 1][i][j] = Xnorm[matlevel][i][j];
				}
			}
		}
		else
		{
			// remove translation and normalize
			matrix[0][3] = 0;
			matrix[1][3] = 0;
			matrix[2][3] = 0;

			scale = sqrt(pow(matrix[0][0], 2.0) + pow(matrix[1][0], 2.0) + pow(matrix[2][0], 2.0));

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					matrix[i][j] = matrix[i][j] / scale;
				}
			}
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					Xnorm[matlevel + 1][i][j] =
						Xnorm[matlevel][i][0] * matrix[0][j]
						+ Xnorm[matlevel][i][1] * matrix[1][j]
						+ Xnorm[matlevel][i][2] * matrix[2][j]
						+ Xnorm[matlevel][i][3] * matrix[3][j];
				}
			}
		}

		matlevel += 1;
		return GZ_SUCCESS;
	}
	else
	{
		return GZ_FAILURE;
	}
}

int GzRender::GzPopMatrix()
{
	/*
	- pop a matrix off the Ximage stack
	- check for stack underflow
	*/
	if (matlevel == 0)
	{
		return GZ_FAILURE;
	}
	else
	{
		matlevel -= 1;
	}

	return GZ_SUCCESS;
}

int GzRender::GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* write pixel values into the buffer */
	// first check if pixel coordinates are valid
	if (i >= 0 && i <= xres - 1 && j >= 0 && j <= yres - 1)
	{
		// clamp RGB values if necessary
		if (r < 0) r = 0;
		if (g < 0) g = 0;
		if (b < 0) b = 0;
		if (r > 4095) r = 4095;
		if (g > 4095) g = 4095;
		if (b > 4095) b = 4095;

		pixelbuffer[ARRAY(i, j)].red = r;
		pixelbuffer[ARRAY(i, j)].green = g;
		pixelbuffer[ARRAY(i, j)].blue = b;
		pixelbuffer[ARRAY(i, j)].alpha = a;
		pixelbuffer[ARRAY(i, j)].z = z;

		return GZ_SUCCESS;
	}
	else return GZ_FAILURE;

}

int GzRender::GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* retrieve a pixel information from the pixel buffer */

	*r = pixelbuffer[ARRAY(i, j)].red;
	*g = pixelbuffer[ARRAY(i, j)].green;
	*b = pixelbuffer[ARRAY(i, j)].blue;
	*a = pixelbuffer[ARRAY(i, j)].alpha;
	*z = pixelbuffer[ARRAY(i, j)].z;

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2File(FILE* outfile)
{
	/* write image to ppm file -- "P6 %d %d 255\r" */

	char rgb_buffer[3];

	// write header
	fprintf(outfile, "P6 %d %d 255\r", xres, yres);

	for (int i = 0; i < xres * yres; i++)
	{
		// write RGB data
		rgb_buffer[0] = pixelbuffer[i].red >> 4;
		rgb_buffer[1] = pixelbuffer[i].green >> 4;
		rgb_buffer[2] = pixelbuffer[i].blue >> 4;

		fwrite(rgb_buffer, sizeof(char), sizeof(rgb_buffer), outfile);
	}

	return GZ_SUCCESS;
}

int GzRender::GzFlushDisplay2FrameBuffer()
{
	/* write pixels to framebuffer:
	- put the pixels into the frame buffer
	- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
	- NOT red, green, and blue !!!
	*/

	for (int i = 0; i < xres * yres * 3; i += 3)
	{
		// write BGR data
		framebuffer[i] = pixelbuffer[i / 3].blue >> 4;
		framebuffer[i + 1] = pixelbuffer[i / 3].green >> 4;
		framebuffer[i + 2] = pixelbuffer[i / 3].red >> 4;
	}

	return GZ_SUCCESS;
}

// helper function to apply transformation to model coordinates
void GzRender::ApplyTransformation(GzCoord v, GzCoord trans_v, GzMatrix stack, bool is_position)
{
	float homo_v[4];
	float xhomo_v[4];
	int cull_count = 0;

	homo_v[X] = v[X];
	homo_v[Y] = v[Y];
	homo_v[Z] = v[Z];
	homo_v[W] = 1;

	// apply Ximage
	xhomo_v[X] = stack[X][X] * homo_v[X] + stack[X][Y] * homo_v[Y] +
		stack[X][Z] * homo_v[Z] + stack[X][W] * homo_v[W];
	xhomo_v[Y] = stack[Y][X] * homo_v[X] + stack[Y][Y] * homo_v[Y] +
		stack[Y][Z] * homo_v[Z] + stack[Y][W] * homo_v[W];
	xhomo_v[Z] = stack[Z][X] * homo_v[X] + stack[Z][Y] * homo_v[Y] +
		stack[Z][Z] * homo_v[Z] + stack[Z][W] * homo_v[W];
	xhomo_v[W] = stack[W][X] * homo_v[X] + stack[W][Y] * homo_v[Y] +
		stack[W][Z] * homo_v[Z] + stack[W][W] * homo_v[W];

	trans_v[X] = xhomo_v[X] / xhomo_v[W];
	trans_v[Y] = xhomo_v[Y] / xhomo_v[W];
	trans_v[Z] = xhomo_v[Z] / xhomo_v[W];

	if (is_position && xhomo_v[Z] < 0)
	{
		cull = true;
	}
}

void GzRender::CalculateC(GzCoord norm, GzColor c)
{
	float N_dot_L, N_dot_E, R_dot_E, scale;
	GzColor s_sum, d_sum;
	GzCoord E = { 0, 0, -1 }, R, N;

	ApplyTransformation(norm, N, Xnorm[matlevel], false);

	for (int i = 0; i < 3; i++)
	{
		s_sum[i] = 0.0;
		d_sum[i] = 0.0;
	}

	for (int i = 0; i < numlights; i++)
	{
		// check signs of N dot L and N dot E
		N_dot_L = dot(N, lights[i].direction);
		N_dot_E = dot(N, E);

		// different signs, skip
		if ((N_dot_L < 0 && N_dot_E > 0) || (N_dot_L > 0 && N_dot_E < 0))
		{
			continue;
		}
		// both negative, flip N
		if (N_dot_L < 0 && N_dot_E < 0)
		{
			N[X] = -(N[X]);
			N[Y] = -(N[Y]);
			N[Z] = -(N[Z]);
		}
		N_dot_L = dot(N, lights[i].direction);

		R[X] = (2 * N_dot_L * N[X]) - lights[i].direction[X];
		R[Y] = (2 * N_dot_L * N[Y]) - lights[i].direction[Y];
		R[Z] = (2 * N_dot_L * N[Z]) - lights[i].direction[Z];
		R_dot_E = dot(R, E);

		// R = 0 to 1
		if (R_dot_E < 0.0)
			R_dot_E = 0.0;
		else if (R_dot_E > 1.0)
			R_dot_E = 1.0;

		for (int j = 0; j < 3; j++)
		{
			s_sum[j] += lights[i].color[j] * pow(R_dot_E, spec);
			d_sum[j] += lights[i].color[j] * N_dot_L;
		}
	}

	for (int i = 0; i < 3; i++)
	{
		// compute color using the shading equation
		c[i] = Ks[i] * s_sum[i] + Kd[i] * d_sum[i] + Ka[i] * ambientlight.color[i];
	}

	return;
}

// helper function to sort vertices in screen coordinates
int* sortVerts(GzPointer *verts)
{
	// sort vertices by their Y values in clockwise direction
	static int ordered_v[3];
	GzCoord *v;
	v = (GzCoord*)verts[0];

	int max_i = 0, mid_i = 0, min_i = 0;
	float max_val = INT_MIN, mid_val = INT_MIN, min_val = INT_MIN;
	float A, B, C, X_val, Y_val, X_comp;

	// sort Y values first
	for (int i = 0; i < 3; i++)
	{
		if (v[i][Y] > max_val)
		{
			min_val = mid_val;
			mid_val = max_val;
			max_val = v[i][Y];
			min_i = mid_i;
			mid_i = max_i;
			max_i = i;
		}
		else if (v[i][Y] > mid_val)
		{
			min_val = mid_val;
			mid_val = v[i][Y];
			min_i = mid_i;
			mid_i = i;
		}
		else
		{
			min_val = v[i][Y];
			min_i = i;
		}
	}

	// cases when 2 verts have the same y value
	if (max_val == mid_val)
	{
		if (v[max_i][X] > v[mid_i][X])
		{
			ordered_v[1] = max_i;
			ordered_v[2] = mid_i;
		}
		else
		{
			ordered_v[1] = mid_i;
			ordered_v[2] = max_i;
		}
		ordered_v[0] = min_i;
	}
	else if (mid_val == min_val)
	{
		if (v[mid_i][X] > v[min_i][X])
		{
			ordered_v[1] = mid_i;
			ordered_v[0] = min_i;
		}
		else
		{
			ordered_v[1] = min_i;
			ordered_v[0] = mid_i;
		}
		ordered_v[2] = max_i;
	}
	else
	{
		/* use edge equation to compute x value of max_vert
		/* when y of max_vert = y_of mid_vert */
		X_val = v[max_i][X];
		Y_val = v[max_i][Y];
		A = v[min_i][Y] - v[max_i][Y];
		B = v[max_i][X] - v[min_i][X];
		C = ((v[min_i][X] - v[max_i][X]) * Y_val) -
			((v[min_i][Y] - v[max_i][Y]) * X_val);

		X_comp = -(B * v[mid_i][Y]) - C;
		X_comp = X_comp / A;

		// vert with the larger x value = 2nd vert counting clockwise
		if (X_comp > v[mid_i][X])
		{
			ordered_v[1] = max_i;
			ordered_v[2] = mid_i;
		}
		else
		{
			ordered_v[1] = mid_i;
			ordered_v[2] = max_i;
		}
		ordered_v[0] = min_i;
	}

	return ordered_v;
}

// helper function to create edge DDA
void createEdges(edge_dda *edge, int v_start_index, int v_end_index,
	GzCoord v[], GzCoord n[], GzColor c[])
{
	/* create edge dda data structure */
	for (int i = 0; i < 3; i++)
	{
		edge->start_v[i] = v[v_start_index][i];
		edge->end_v[i] = v[v_end_index][i];
		edge->current_v[i] = v[v_start_index][i];

		edge->start_n[i] = n[v_start_index][i];
		edge->end_n[i] = n[v_end_index][i];
		edge->current_n[i] = n[v_start_index][i];

		edge->start_c[i] = c[v_start_index][i];
		edge->end_c[i] = c[v_end_index][i];
		edge->current_c[i] =c[v_start_index][i];
	}

	edge->slope_x = (edge->end_v[X] - edge->start_v[X]) /
		(edge->end_v[Y] - edge->start_v[Y]);
	edge->slope_z = (edge->end_v[Z] - edge->start_v[Z]) /
		(edge->end_v[Y] - edge->start_v[Y]);

	for (int i = 0; i < 3; i++)
	{
		edge->slope_n[i] = (edge->end_n[i] - edge->start_n[i]) /
			(edge->end_v[Y] - edge->start_v[Y]);
		edge->slope_c[i] = (edge->end_c[i] - edge->start_c[i]) /
			(edge->end_v[Y] - edge->start_v[Y]);
	}

	return;
}

// helper function to update edge DDA
void updateEdges(edge_dda *edge, float delta_y)
{
	/* modify values in edge dda*/
	edge->current_v[X] += ((edge->slope_x) * delta_y);
	edge->current_v[Y] += delta_y;
	edge->current_v[Z] += ((edge->slope_z) * delta_y);

	for (int i = 0; i < 3; i++)
	{
		edge->current_n[i] += ((edge->slope_n[i]) * delta_y);
		edge->current_c[i] += ((edge->slope_c[i]) * delta_y);
	}

	return;
}

// helper function to create span line DDA
void createSpan(span_dda *span, edge_dda *left, edge_dda *right)
{
	/* create span line dda data structure */
	span->left_x = left->current_v[X];
	span->left_z = left->current_v[Z];
	span->right_x = right->current_v[X];
	span->right_z = right->current_v[Z];
	span->current_x = span->left_x;
	span->current_z = span->left_z;
	span->slope_z = (span->right_z - span->left_z) / (span->right_x - span->left_x);

	for (int i = 0; i < 3; i++)
	{
		span->left_n[i] = left->current_n[i];
		span->right_n[i] = right->current_n[i];
		span->current_n[i] = span->left_n[i];
		span->slope_n[i] = (span->right_n[i] - span->left_n[i]) / (span->right_x - span->left_x);

		span->left_c[i] = left->current_c[i];
		span->right_c[i] = right->current_c[i];
		span->current_c[i] = span->left_c[i];
		span->slope_c[i] = (span->right_c[i] - span->left_c[i]) / (span->right_x - span->left_x);
	}

	return;
}

// helper function to update span line DDA
void updateSpan(span_dda *span, float delta_x)
{
	/* modify values in span line dda*/
	span->current_x += delta_x;
	span->current_z += ((span->slope_z) * delta_x);
	for (int i = 0; i < 3; i++)
	{
		span->current_n[i] += ((span->slope_n[i]) * delta_x);
		span->current_c[i] += ((span->slope_c[i]) * delta_x);
	}

	return;
}

int GzRender::GzPutAttribute(int numAttributes, GzToken	*nameList, GzPointer *valueList)
{
	/*
	-- Set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
	-- In later homeworks set shaders, interpolaters, texture maps, and lights
	*/

	for (int i = 0; i < numAttributes; i++)
	{
		if (nameList[i] == GZ_RGB_COLOR)
		{
			float *color;
			color = (float*)valueList[i];

			for (int i = 0; i < 3; i++)
			{
				flatcolor[i] = color[i];
			}
		}
		else if (nameList[i] == GZ_INTERPOLATE)
		{
			int *interp;
			interp = (int*)valueList[i];

			interp_mode = *interp;

		}
		else if (nameList[i] == GZ_DIRECTIONAL_LIGHT)
		{
			float *d_light;
			d_light = (float*)valueList[i];

			if (numlights >= 0 && numlights < 10)
			{
				lights[numlights].direction[0] = d_light[0];
				lights[numlights].direction[1] = d_light[1];
				lights[numlights].direction[2] = d_light[2];
				lights[numlights].color[0] = d_light[3];
				lights[numlights].color[1] = d_light[4];
				lights[numlights].color[2] = d_light[5];

				numlights++;
			}
		}
		else if (nameList[i] == GZ_AMBIENT_LIGHT)
		{
			float *a_light;
			a_light = (float*)valueList[i];

			ambientlight.direction[0] = a_light[0];
			ambientlight.direction[1] = a_light[1];
			ambientlight.direction[2] = a_light[2];
			ambientlight.color[0] = a_light[3];
			ambientlight.color[1] = a_light[4];
			ambientlight.color[2] = a_light[5];

		}
		else if (nameList[i] == GZ_AMBIENT_COEFFICIENT)
		{
			float *a_c;
			a_c = (float*)valueList[i];
			for (int i = 0; i < 3; i++)
			{
				Ka[i] = a_c[i];
			}
		}
		else if (nameList[i] == GZ_DIFFUSE_COEFFICIENT)
		{
			float *d_c;
			d_c = (float*)valueList[i];
			for (int i = 0; i < 3; i++)
			{
				Kd[i] = d_c[i];
			}
		}
		else if (nameList[i] == GZ_SPECULAR_COEFFICIENT)
		{
			float *s_c;
			s_c = (float*)valueList[i];
			for (int i = 0; i < 3; i++)
			{
				Ks[i] = s_c[i];
			}
		}
		else if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT)
		{
			float *s;
			s = (float*)valueList[i];

			spec = *s;
		}
		else
		{
			return GZ_FAILURE;
		}
	}
	return GZ_SUCCESS;
}

int GzRender::GzPutTriangle(int	numParts, GzToken *nameList, GzPointer *valueList)
/* numParts - how many names and values */
{
	/*
	-- Pass in a triangle description with tokens and values corresponding to
	GZ_NULL_TOKEN:		do nothing - no values
	GZ_POSITION:		3 vert positions in model space
	-- Invoke the rastrizer/scanline framework
	-- Return error code
	*/

	GzCoord transformed_v[3], interp_n;
	GzCoord n[3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	GzColor c[3] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	GzColor color;
	GzCoord *ptr;
	GzPointer v[3];

	int *sorted_v;
	edge_dda e12, e23, e13;
	span_dda span_line;
	edge_dda *e_12 = &e12;
	edge_dda *e_23 = &e23;
	edge_dda *e_13 = &e13;
	edge_dda *e_long, *e_short, *e_left, *e_right;
	span_dda *span = &span_line;

	int idx = 0, v_1_index, v_2_index, v_3_index, new_Z;
	float scale;

	for (int i = 0; i < numParts; i++)
	{
		ptr = (GzCoord*)valueList[i];

		if (nameList[i] == GZ_POSITION)
		{
			// apply transformation to convert model coordinates to screen
			for (int i = 0; i < 3; i++)
			{
				ApplyTransformation(ptr[i], transformed_v[i], Ximage[matlevel], true);
			}
		}
		else if (nameList[i] == GZ_NORMAL)
		{
			for (int i = 0; i < 3; i++)
			{
				n[i][X] = ptr[i][X];
				n[i][Y] = ptr[i][Y];
				n[i][Z] = ptr[i][Z];
			}
		}
	}

	if (interp_mode == GZ_FLAT || interp_mode == GZ_COLOR)
	{
		for (int i = 0; i < 3; i++)
		{
			// calculate color at each vertex
			CalculateC(n[i], c[i]);
		}
		// use first normal for flat color
		for (int i = 0; i < 3; i++)
		{
			flatcolor[i] = c[0][i];
		}
	}

	// check for Z < 0
	if (cull == false)
	{
		// sort
		v[0] = (GzPointer)transformed_v;
		sorted_v = sortVerts(v);
		v_1_index = *sorted_v;
		v_2_index = *(sorted_v + 1);
		v_3_index = *(sorted_v + 2);

		// make edge DDAs
		createEdges(e_12, v_1_index, v_2_index, transformed_v, n, c);
		createEdges(e_13, v_1_index, v_3_index, transformed_v, n, c);

		// keep track of long, short, left, and right edge for traversal
		e_left = &e13;
		e_right = &e12;

		if (e13.end_v[Y] > e12.end_v[Y])
		{
			e_long = &e13;
			e_short = &e12;
			createEdges(e_23, v_2_index, v_3_index, transformed_v, n, c);
		}
		else
		{
			e_long = &e12;
			e_short = &e13;
			createEdges(e_23, v_3_index, v_2_index, transformed_v, n, c);
		}

		// traverse down to nearest integer y point
		updateEdges(e_long, (ceil(e_long->start_v[Y]) -
			e_long->start_v[Y]));
		updateEdges(e_short, (ceil(e_short->start_v[Y]) -
			e_short->start_v[Y]));

		// traverse to end of long edge / border
		while (e_long->current_v[Y] < e_long->end_v[Y] &&
			e_long->current_v[Y] <= (yres - 1))
		{
			// switch to 2nd short edge when end of 1st short edge is reached
			if (e_short->current_v[Y] > e_short->end_v[Y])
			{
				e_short = &e23;

				if (e_long == &e13)
				{
					e_right = &e23;
				}
				else e_left = &e23;

				updateEdges(e_short, (ceil(e_short->start_v[Y]) -
					e_short->start_v[Y]));
			}

			// make span line
			createSpan(span, e_left, e_right);

			// traverse right to nearnest interger x point
			updateSpan(span, (ceil(span_line.left_x) - span_line.left_x));

			// traverse to end of span line / border
			while (span_line.current_x < span_line.right_x &&
				span_line.current_x <= (xres - 1))
			{
				if (e_long->current_v[Y] >= 0 && span_line.current_x >= 0
					&& e_long->current_v[Y] <= (yres - 1)
					&& span_line.current_x <= (xres - 1))
				{
					// flat shading
					if (interp_mode == GZ_FLAT)
					{
						color[RED] = ctoi(flatcolor[RED]);
						color[GREEN] = ctoi(flatcolor[GREEN]);
						color[BLUE] = ctoi(flatcolor[BLUE]);
					}
					// Gouraud shading
					else if (interp_mode == GZ_COLOR)
					{
						color[RED] = ctoi(span_line.current_c[RED]);
						color[GREEN] = ctoi(span_line.current_c[GREEN]);
						color[BLUE] = ctoi(span_line.current_c[BLUE]);
					}
					//Phong shading
					else if (interp_mode == GZ_NORMALS)
					{
						// normalize N first
						for (int i = 0; i < 3; i++)
						{
							interp_n[i] = span_line.current_n[i];
						}

						scale = sqrt(pow(interp_n[X], 2.0) +
							pow(interp_n[Y], 2.0) + pow(interp_n[Z], 2.0));

						for (int i = 0; i < 3; i++)
						{
							interp_n[i] = interp_n[i] / scale;
						}

						// calculate color for each pixel
						CalculateC(interp_n, color);
						color[RED] = ctoi(color[RED]);
						color[GREEN] = ctoi(color[GREEN]);
						color[BLUE] = ctoi(color[BLUE]);
					}

					new_Z = span_line.current_z;

					// only write pixel if Z is smaller (closer)
					if (new_Z < pixelbuffer[ARRAY((int)span_line.current_x,
						(int)e_long->current_v[Y])].z)
					{
						GzPut(span_line.current_x, e_long->current_v[Y],
							color[RED], color[GREEN], color[BLUE], 1, new_Z);
					}
				}
				updateSpan(span, 1.0);
			}
			updateEdges(e_long, 1.0);
			updateEdges(e_short, 1.0);
		}
	}
	else
	{
		cull = false;
	}

	return GZ_SUCCESS;
}
