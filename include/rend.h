#include	"gz.h"
#ifndef GZRENDER_
#define GZRENDER_

/* Camera defaults */
#define	DEFAULT_FOV		35.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	  {0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	  {0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		  32

#define	MATLEVELS	  100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#ifndef GzTexture
#define GzTexture	GzPointer
#endif

#define W	3

struct edge_dda {
	GzCoord start_v;
	GzCoord end_v;
	GzCoord current_v;
	GzCoord start_n;
	GzCoord end_n;
	GzCoord current_n;
	GzColor start_c;
	GzColor end_c;
	GzColor current_c;
	float slope_x;
	float slope_z;
	GzCoord slope_n;
	GzColor slope_c;

};

struct span_dda {
	float left_x;
	float left_z;
	float right_x;
	float right_z;
	float current_x;
	float current_z;
	float slope_z;
	GzCoord left_n;
	GzCoord right_n;
	GzCoord current_n;
	GzCoord slope_n;
	GzColor left_c;
	GzColor right_c;
	GzColor current_c;
	GzColor slope_c;
};

class GzRender {			/* define a renderer */

  public:
  	unsigned short	xres;
  	unsigned short	yres;
  	GzPixel		      *pixelbuffer;		/* frame buffer array */
  	char*           framebuffer;

  	GzCamera		m_camera;
  	short		    matlevel;	        /* top of stack - current xform */
  	GzMatrix		Ximage[MATLEVELS];	/* stack of xforms (Xsm) */
  	GzMatrix		Xnorm[MATLEVELS];	/* xforms for norms (Xim) */
  	GzMatrix		Xsp;		        /* NDC to screen (pers-to-screen) */
  	GzColor		  flatcolor;          /* color state for flat shaded triangles */
  	int         interp_mode;
  	int			    numlights;
  	GzLight		  lights[MAX_LIGHTS];
  	GzLight		  ambientlight;
  	GzColor		  Ka, Kd, Ks;
  	float		    spec;		/* specular power */
  	GzTexture		tex_fun;    /* tex_fun(float u, float v, GzColor color) */


    // Constructors
  	GzRender(int xRes, int yRes);
  	~GzRender();

  	// Display methods
  	int GzDefault();
  	int GzBeginRender();
  	int GzPut(int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z);
  	int GzGet(int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth	*z);

  	int GzFlushDisplay2File(FILE* outfile);
  	int GzFlushDisplay2FrameBuffer();

  	// Render methods
  	int GzPutAttribute(int numAttributes, GzToken *nameList, GzPointer *valueList);
  	int GzPutTriangle(int numParts, GzToken *nameList, GzPointer *valueList);

  	//
  	int GzDefaultCamera();
  	int GzPutCamera(GzCamera camera);
  	int GzPushMatrix(GzMatrix matrix);
  	int GzPopMatrix();
  	void ApplyTransformation(GzCoord v, GzCoord trans_v, GzMatrix stack, bool is_position);
  	void CalculateC(GzCoord norm, GzColor c);

  	// Extra methods: NOT part of API - just for general assistance */
  	inline int ARRAY(int x, int y) {return (x+y*xres);}	/* simplify fbuf indexing */
  	inline short	ctoi(float color) {return(short)((int)(color * ((1 << 12) - 1)));}		/* convert float color to GzIntensity short */

  	// Object Translation
  	int GzRotXMat(float degree, GzMatrix mat);
  	int GzRotYMat(float degree, GzMatrix mat);
  	int GzRotZMat(float degree, GzMatrix mat);
  	int GzTrxMat(GzCoord translate, GzMatrix mat);
  	int GzScaleMat(GzCoord scale, GzMatrix mat);

};

#endif
