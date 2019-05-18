#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <windows.h>
#endif

//#include <GL/gl.h>
//#include <GL/glut.h>
#include <GLUT/glut.h> // OS X

#ifndef M_PI
#define M_PI 3.141592653589793238462643
#endif

#define min(a, b) ((a)<(b)?(a):(b))
#define max(a, b) ((a)>(b)?(a):(b))
#define clamp(x, a, b) min(max(x, a), b)

#define WINDOW_TITLE "I3D Assignment 3"

enum RenderOptions {
	OPT_WIREFRAME,
	OPT_RENDER_SATURN,
	OPT_RENDER_RINGS,
	OPT_DRAW_AXES,
	OPT_ORTHOGRAPHIC,
	
	OPT_ANIMATION,
	OPT_LIGHTING,
	OPT_CULLING,
	OPT_CULL_FACE,
	OPT_DRAW_NORMALS,
	
	OPTSIZE //NOTE: this must be the last enumerator.
};

typedef struct Vec3f {
	float x, y, z;
} Vec3f;

typedef struct TexUV {
	float u, v;
} TexUV;

typedef struct Camera {
	Vec3f rot; // Euler angles - pitch/yaw/roll.
	Vec3f pos;
	float zoom;
	float sensitivity;
	float fov;
	bool rotating;
	bool zooming;
	bool panning;
} Camera;

typedef struct Mesh {
	Vec3f** grid;
	Vec3f** normals;
	int rows;
	int cols;
	int type;
} Mesh;

static int primitiveMode = 0;
static const int tessSaturn = 16;
static const int tessParticle = 8;
static int tessDisc = 4;
static const int maxTes = 1024;
static const int minTes = 4;
static bool options[OPTSIZE];
static Camera camera;
static Mesh saturn;
static Mesh rings;
static Mesh discParticle;
static float saturnTilt = 0.0f;
static const float clipFar = (3560000 * 2.1);
static const float clipNear = (3560000 * 2.1) / 10000.0;

#define SATURN_DAY (((10*60*60) + (32*60) + (35)) * 1000)

static float saturn_rotation = 0;
static float saturn_rotate_velocity = 0.012;

// lighting 
static float light0_ambient[] = {0.1, 0.1, 0.1, 1.0};
static float light0_diffuse[] = {1.0, 1.0, 1.0, 1.0};
static float light0_specular[] = {1.0, 1.0, 1.0, 1.0};
static float light0_position[] = {1.0, 1.0, 1.0, 0.0};
static const float shininess = {50};
const float white[] = {1.0, 1.0, 1.0};
const float blue[] = {0.0, 0.0, 1.0};
// *end

// osd
/* set update interval to 1 sec */
static const int osdUpdateInterval = 1000;
static char fps[50];
static char frameTime[50];
static char particles[50];
static int frames;
static int particleCount;
static float osdElapsedTime;
// *end


#define SATURN 1
#define DISC 2
#define PARTICLE 3
#define SATURN_RADIUS 60268
#define RING_DIST_INNER 74510
#define RING_DIST_OUTER 140390
#define SATURN_TILT_MAX 27
#define SATURN_SCALE_HEIGHT 0.9 // Squash saturn by 10%.

//#define DEBUG_MESH
//#define DEBUG_CAMERA

static const GLenum primitiveTypes[] = {GL_QUADS, GL_POINTS};
static bool primitiveTypeSphere = false;

#define checkGLErrors() _checkGLErrors(__LINE__)
void _checkGLErrors(int line)
{
	GLuint err;
	if ((err = glGetError()) != GL_NO_ERROR)
		fprintf(stderr, "OpenGL error caught on line %i: %s\n", line, gluErrorString(err));
}

void drawAxes(float len, int lineWidth)
{
	if (!options[OPT_DRAW_AXES])
		return;
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT);
	glLineWidth(lineWidth);
	glBegin(GL_LINES);
	glColor3f(1, 0, 0); glVertex3f(0, 0, 0); glVertex3f(len, 0, 0);
	glColor3f(0, 1, 0); glVertex3f(0, 0, 0); glVertex3f(0, len, 0);
	glColor3f(0, 0, 1); glVertex3f(0, 0, 0); glVertex3f(0, 0, len);
	glEnd();
	glPopAttrib();
}

void allocateMesh(Mesh* mesh)
{
	// Allocate grid.
	mesh->grid = (Vec3f**)malloc(sizeof(Vec3f*) * mesh->rows); // Allocate array of pointers.
	for (int i = 0; i < mesh->rows; i++)
		mesh->grid[i] = (Vec3f*)malloc(sizeof(Vec3f) * mesh->cols); // Allocate arrays of data.

	// Allocate Normal grid.
	mesh->normals = (Vec3f**)malloc(sizeof(Vec3f*) * mesh->rows); // Allocate array of pointers.
	for (int i = 0; i < mesh->rows; i++)
		mesh->normals[i] = (Vec3f*)malloc(sizeof(Vec3f) * mesh->cols); // Allocate arrays of normal data.
}

void createSphere(Mesh *mesh, int cols, int rows, float r, int type)
{
#ifdef DEBUG_MESH
	fprintf(stderr, "new sphere %i %i\n", cols, rows);
#endif
	// Allocate a mesh to return.
	mesh->rows = rows;
	mesh->cols = cols;
	mesh->type = type; 

	allocateMesh(mesh);
	
	// Create the vertices for a sphere.
	for (int i = 0; i < mesh->rows; ++i) 
	{
		float v = M_PI * (float)i / (mesh->rows-1);
		for (int j = 0; j < mesh->cols; ++j)
		{
			float u = 2.0 * M_PI * (float)j / (mesh->cols-1);
			
			// Parametric coordinates.
			//http://en.wikipedia.org/wiki/Parametric_surface
			mesh->grid[i][j].x = r * cos(u) * sin(v);
			mesh->grid[i][j].y = r * sin(u) * sin(v);
			mesh->grid[i][j].z = r * cos(v);
			
			mesh->normals[i][j].x = (mesh->grid[i][j].x)/r;
			mesh->normals[i][j].y = (mesh->grid[i][j].y)/r;
			mesh->normals[i][j].z = (mesh->grid[i][j].z)/r;
		}
	}
}

void createDisc(Mesh *mesh, int cols, int rows, float r, float R, int type)
{
#ifdef DEBUG_MESH
	fprintf(stderr, "new disc %i %i\n", cols, rows);
#endif
	// Allocate a mesh to return.
	mesh->rows = rows;
	mesh->cols = cols;
	mesh->type = type; 

	allocateMesh(mesh);	
	
	// Create the vertices for a disc.
	float d = R - r;
	
	for (int i = 0; i < mesh->rows; ++i)
	{
		float v = (float)i / (mesh->rows-1);
		for (int j = 0; j < mesh->cols; ++j)
		{
			float u = 2.0 * M_PI * (float)j / (mesh->cols-1);
			mesh->grid[i][j].x = cos(u) * (r + d * v);
			mesh->grid[i][j].y = sin(u) * (r + d * v);
			mesh->grid[i][j].z = 0.0;

			mesh->normals[i][j].x = (mesh->grid[i][j].x)/(r + d * v);
			mesh->normals[i][j].y = (mesh->grid[i][j].y)/(r + d * v);
			mesh->normals[i][j].z = (mesh->grid[i][j].z)/(r + d * v);
		}
	}
}

void drawNormalMesh(Mesh* mesh)
{
	if (!options[OPT_DRAW_NORMALS])
		return;
	glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT);
	glColor3fv(blue);
	glBegin(GL_LINES);
	for (int i = 0; i < mesh->rows; i++) {
		for (int j = 0; j < mesh->cols; j++) {
			if (mesh->type == SATURN) {
			glVertex3f(mesh->normals[i][j].x * SATURN_RADIUS, mesh->normals[i][j].y * SATURN_RADIUS,
			mesh->normals[i][j].z * SATURN_RADIUS);  
			glVertex3f(mesh->normals[i][j].x * 1.3 * SATURN_RADIUS , mesh->normals[i][j].y * 1.3
			* SATURN_RADIUS, mesh->normals[i][j].z * 1.3 * SATURN_RADIUS);
			}
			else {
			glVertex3f(mesh->normals[i][j].x * RING_DIST_OUTER, mesh->normals[i][j].y * RING_DIST_OUTER,
			mesh->normals[i][j].z * RING_DIST_OUTER);  
			glVertex3f(mesh->normals[i][j].x * 1.3 * RING_DIST_OUTER, mesh->normals[i][j].y * 1.3
			* RING_DIST_OUTER, mesh->normals[i][j].z * 1.3 * RING_DIST_OUTER);
			}
		}
	}
	glEnd();
	glPopAttrib();
}

void drawMesh(Mesh* mesh)
{
	if (options[OPT_DRAW_NORMALS])
		drawNormalMesh(mesh);

	// Draws each row of the mesh using the set primitive type.
	GLuint primitiveType = primitiveTypes[primitiveMode];
	bool strips = (primitiveType == GL_QUAD_STRIP || primitiveType == GL_TRIANGLE_STRIP);
	bool triangles = (primitiveType == GL_TRIANGLES);
	bool points = (primitiveType == GL_POINTS);

	if (mesh->type == SATURN || mesh->type == PARTICLE)
	{
		primitiveType = GL_QUAD_STRIP;
		strips = true; triangles = false; points = false;
	}

	if (points)
		glPointSize(2.0f);
	else
		glPointSize(1.0f);

	for (int i = 0; i < mesh->rows - 1; ++i)
	{
		glBegin(primitiveType);
		for (int j = 0; j < mesh->cols - 1; ++j)
		{
			const Vec3f a = mesh->grid[i][j];
			const Vec3f b = mesh->grid[i+1][j];
			const Vec3f an = mesh->normals[i][j];
			const Vec3f bn = mesh->normals[i+1][j];

			glNormal3f(an.x, an.y, an.z);
			glVertex3f(a.x, a.y, a.z); // Always draw the first two vertices at this iteration.
			glNormal3f(bn.x, bn.y, bn.z);
			glVertex3f(b.x, b.y, b.z);
			if (!strips)
			{
				// If using quads, another two vertices are needed.
				const Vec3f c = mesh->grid[i][j+1];
				const Vec3f d = mesh->grid[i+1][j+1];
				const Vec3f cn = mesh->normals[i][j+1];
				const Vec3f dn = mesh->normals[i+1][j+1];

				glNormal3f(dn.x, dn.y, dn.z);
				glVertex3f(d.x, d.y, d.z);
				if (triangles)
				{
					// If using triangles, two are needed. The first is already drawn so fill in the second.
					glNormal3f(an.x, an.y, an.z);
					glVertex3f(a.x, a.y, a.z);
					glNormal3f(dn.x, dn.y, dn.z);
					glVertex3f(d.x, d.y, d.z);
				}
				glNormal3f(cn.x, cn.y, cn.z);
				glVertex3f(c.x, c.y, c.z);
			}
		}
		if (strips)
		{
			const Vec3f a = mesh->grid[i][mesh->cols - 1];
			const Vec3f b = mesh->grid[i+1][mesh->cols - 1];
			const Vec3f an = mesh->normals[i][mesh->cols - 1];
			const Vec3f bn = mesh->normals[i+1][mesh->cols - 1];

			glNormal3f(an.x, an.y, an.z);
			glVertex3f(a.x, a.y, a.z); // If using strips, the last two vertices need to be drawn for this stack.
			glNormal3f(bn.x, bn.y, bn.z);
			glVertex3f(b.x, b.y, b.z);
		}
		glEnd();
	}
}

void drawSphereParticleDisc(Mesh* discMesh, Mesh* particleMesh)
{
	for (int i = 0; i < discMesh->rows; ++i)
	{
		for (int j = 0; j < discMesh->cols; ++j)
		{
			glPushMatrix();
			glTranslatef(discMesh->grid[i][j].x, discMesh->grid[i][j].y, discMesh->grid[i][j].z);
			drawMesh(particleMesh);
			glPopMatrix();
		}
	}
}

void drawEPM(Mesh* mesh)
{
	glPushMatrix();
	glColor3f(1, 1, 0);
	glLineWidth(2);

	// Equator
	glBegin(GL_LINE_STRIP);
	for (int i = 0; i < mesh->rows; ++i)
	{
		const int j = mesh->cols/2;
		const Vec3f a = mesh->grid[i][j];
		glVertex3f(a.x, a.y, a.z);
	}
	glEnd();
	
	// Prime Meridian
	glBegin(GL_LINE_STRIP);
	for (int j = 0; j < mesh->cols; ++j)
	{
		const int i = mesh->rows/2;
		const Vec3f a = mesh->grid[i][j];
		glVertex3f(a.x, a.y, a.z);
	}
	glEnd();

	glPopMatrix();
	glLineWidth(1);
}

void freeMesh(Mesh* mesh)
{
	// Free allocated data.
	if (mesh->grid)
	{
		for (int i = 0; i < mesh->rows; ++i)
			free(mesh->grid[i]);
		free(mesh->grid);
	}
	if (mesh->normals)
	{
		for (int i = 0; i < mesh->rows; ++i)
			free(mesh->normals[i]);
		free(mesh->normals);
	}

	mesh->grid = NULL;
	mesh->normals = NULL;
	mesh->rows = 0;
	mesh->cols = 0;
}

void updateRenderState()
{
	// Set wireframe state.
	if (options[OPT_WIREFRAME])
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	else
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void updateLightingAndStuff()
{
	// Set Culling
	if (options[OPT_CULLING])
		glEnable(GL_CULL_FACE);
	else
		glDisable(GL_CULL_FACE);

	// Set Cull Face
	if (options[OPT_CULL_FACE])
		glCullFace(GL_FRONT);
	else
		glCullFace(GL_BACK);

	// Set Lighting
	if (options[OPT_LIGHTING])
		glEnable(GL_LIGHTING);
	else
		glDisable(GL_LIGHTING);
}

void createGeometry()
{
	// Free the meshes (if they exist).
	freeMesh(&saturn);
	freeMesh(&rings);
	freeMesh(&discParticle);

	createSphere(&saturn, tessSaturn + 1, tessSaturn, SATURN_RADIUS, SATURN);
	createDisc(&rings,tessDisc + 1 , tessDisc, RING_DIST_INNER, RING_DIST_OUTER, DISC);
	createSphere(&discParticle, tessParticle + 1, tessParticle, 500.0, PARTICLE);
}

void init()
{
	const float white[] = {1.0, 1.0, 1.0};
	// Initialize options array.
	memset(options, 0, sizeof(bool) * OPTSIZE);
	
	// Initialize camera.
	camera.rot.x = 45.0f; camera.rot.y = 0.0f; camera.rot.z = 0.0f;
	camera.pos.x = 0.0f; camera.pos.y = 0.0f; camera.pos.z = 0.0f;
	camera.zoom = SATURN_RADIUS * 4.0f;
	camera.rotating = camera.zooming = camera.panning = false;
	camera.sensitivity = 0.4;
	camera.fov = 60.0;

	// Initialize geometry.
	createGeometry();

	// Initialize lighting.
	glShadeModel(GL_SMOOTH);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_LIGHT0);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, white);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shininess);

	// osd defaults
	frames = 0;
	particleCount = 0;
	osdElapsedTime = 0;

	// display defaults
	options[OPT_RENDER_SATURN] = true;
	options[OPT_RENDER_RINGS] = true;
}

void reshape(int x, int y)
{
	// Set the portion of the window to draw to.
	glViewport(0, 0, x, y);
	// Modify the projection matrix.
	glMatrixMode(GL_PROJECTION);
	// Clear previous projection.
	glLoadIdentity();
	// Create new projection.
	float aspect = x / (float)y;
	if (options[OPT_ORTHOGRAPHIC])
	{
		float d = clipNear + sqrt(clipFar - clipNear);
		glOrtho(-d * 0.5, d * 0.5, - d * 0.5 / aspect, d * 0.5 / aspect, clipNear, clipFar);
	} 
	else
		gluPerspective(camera.fov, aspect, clipNear, clipFar);
	//IMPORTANT: go back to modifying the modelview matrix.
	glMatrixMode(GL_MODELVIEW);
}

void setupCamera()
{
	// Handle orthographic projection "zooming".
	if (options[OPT_ORTHOGRAPHIC])
	{
		float s = 1000.0/((camera.zoom - clipNear)*tan(camera.fov));
		glTranslatef(0, 0, -clipFar/2.0);
		glScalef(s, s, 1.0);
	} 
	else
		glTranslatef(0, 0, -camera.zoom);


	glTranslatef(-camera.pos.x, -camera.pos.y, -camera.pos.z);

	glRotatef(camera.rot.x, 1, 0, 0);
	glRotatef(camera.rot.y, 0, 1, 0);
	glRotatef(camera.rot.z, 0, 0, 1);
}

void drawEntity()
{
	const float lightBrown[] = {224.0/255.0, 188.0/255.0, 126.0/255.0};
	const float darkBrown[] = {120.0/255.0, 108.0/255.0, 86.0/255.0};

	// Draw saturn.
	glPushMatrix();
	glScalef(1.0, 1.0, SATURN_SCALE_HEIGHT);
		
	glRotatef(saturn_rotation, 0.0, 0.0, 1.0);

	if (!options[OPT_WIREFRAME]) {
		glColor3fv(lightBrown);
	}
	if (options[OPT_RENDER_SATURN])
		drawMesh(&saturn);

	if (!options[OPT_WIREFRAME]) {
		glColor3fv(darkBrown);
	}
	if (options[OPT_RENDER_RINGS]) {
		if (!primitiveTypeSphere)
			drawMesh(&rings);
		else
			drawSphereParticleDisc(&rings, &discParticle);
	}

#if 0
	if (options[OPT_WIREFRAME])
		drawEPM(&saturn); // Draw equator and meridian for Saturn.
#endif

	glPopMatrix();
}

void osd()
{
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
		glLoadIdentity();
		glOrtho(0, glutGet(GLUT_WINDOW_WIDTH)/2, 0, glutGet(GLUT_WINDOW_HEIGHT)/2, 0, 10000);
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();
		glColor3f(1.0, 1.0, 0.0);
		glRasterPos2f(10, 50);
		for (int i = 0; i < strlen(fps); i++)
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, fps[i]);
		glRasterPos2f(10, 35);
		for (int i = 0; i < strlen(frameTime); i++)
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, frameTime[i]);
		glRasterPos2f(10, 20);
		for (int i = 0; i < strlen(particles); i++)
			glutBitmapCharacter(GLUT_BITMAP_9_BY_15, particles[i]);	
		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
}

void display()
{
	glLoadIdentity();
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// light
	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);

	// Position the camera
	setupCamera();

	// update Lighting, culling, etc.
	updateLightingAndStuff();

	glRotatef(-90.0 + saturnTilt, 1.0, 0.0, 0.0); // Tilt of planet.

	if (options[OPT_DRAW_AXES])
		drawAxes(20000.0, 2);

	glColor3f(1, 1, 1);
	
	// Draw saturn and all moons.
	drawEntity();

	// OSD
	osd();

	// Display result.
	glutSwapBuffers();

	// Check for gl errors.
	checkGLErrors();
}

void getParticleCount(int* particleCount)
{
	if (options[OPT_RENDER_RINGS])
		*particleCount += (tessDisc)*(tessDisc);

}

void idle()
{	
	static int now_time, last_time = 0, elapsed_time;

	now_time = glutGet(GLUT_ELAPSED_TIME);
	elapsed_time = now_time - last_time;
	last_time = now_time;

	// Increment frames
	frames++;
	// Update osd elapsed time	
	osdElapsedTime += elapsed_time;
	if (osdElapsedTime >= osdUpdateInterval)
	{
		getParticleCount(&particleCount);
		sprintf(fps, "%-25s%i", "Frame Rate (fps) : ", frames);
		sprintf(frameTime, "%-25s%i", "Frame Time (ms) : ", (int)(osdElapsedTime/(frames)));
		sprintf(particles, "%-25s%i", "Particles : ", particleCount);		
		osdElapsedTime = 0;
		frames = 0;
		particleCount = 0;
	}

	if (options[OPT_ANIMATION]) {
		saturn_rotation += elapsed_time * saturn_rotate_velocity;
		if (saturn_rotation > 360)
			saturn_rotation -= 360;
	}

	glutPostRedisplay();
}

void cleanup()
{
	freeMesh(&saturn);
	freeMesh(&rings);
	freeMesh(&discParticle);
}


void mouseDown(int button, int state, int x, int y)
{
#ifdef DEBUG_CAMERA
	printf("%d %d\n", camera.rotating, camera.zooming);
#endif
	if (button == GLUT_LEFT_BUTTON)
		camera.rotating = state == GLUT_DOWN;
	if (button == GLUT_RIGHT_BUTTON)
		camera.zooming = state == GLUT_DOWN;
	if (button == GLUT_MIDDLE_BUTTON)
		camera.panning = state == GLUT_DOWN;
}

void mouseMove(int x, int y)
{
	static int prevX, prevY;
	
	int movX = x - prevX;
	int movY = y - prevY;
	
	if (camera.rotating)
	{
		// Rotate the camera.
		camera.rot.x += movY * camera.sensitivity;
		camera.rot.y += movX * camera.sensitivity;
		
		// Clamp camera's pitch.
		camera.rot.x = clamp(camera.rot.x, -90.0, 90.0);
	}
	
	if (camera.zooming)
	{
		// Zoom in/out (speed is proportional to the current zoom).
		camera.zoom -= movY * camera.zoom * camera.sensitivity * 0.1;
		
		// Clamp the camera's zoom.
		camera.zoom = clamp(camera.zoom, clipNear, clipFar / 2.0);
	}

	if (camera.panning)
	{
		// pan along x axis (speed proportional to zoom and window size)
		camera.pos.x += (-movX  * camera.zoom * camera.sensitivity) / (glutGet(GLUT_WINDOW_WIDTH)/2);
		// pan along y axis (speed proportional to zoom and window size)
		camera.pos.y += (movY  * camera.zoom * camera.sensitivity) / (glutGet(GLUT_WINDOW_HEIGHT)/2);
	}
	
	prevX = x;
	prevY = y;
}

void keyDown(unsigned char k, int x, int y)
{

	switch (k) {
	case 27:
	case 'q':
		// Exit if esc or q is pressed.
		cleanup();
		exit(0);
		break;
	case 'w':
		options[OPT_WIREFRAME] = !options[OPT_WIREFRAME];
		break;
	case 'a':
		options[OPT_DRAW_AXES] = !options[OPT_DRAW_AXES];
		break;
	case 'p':
		options[OPT_ORTHOGRAPHIC] = !options[OPT_ORTHOGRAPHIC];
		reshape(glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT));
		break;
	case '+':
	case '-':
		if (k == '+')
			tessDisc = min(maxTes, tessDisc * 2.0);
		else
			tessDisc = max(minTes, tessDisc / 2.0);
		createGeometry();
		break;
	case 'i':
	case 'I':
		if (k == 'I')
		{
			if (saturnTilt < SATURN_TILT_MAX)
				saturnTilt += 1.0;
		}
		else 
		{
			saturnTilt -= 1.0;
			saturnTilt = fmin(fmax(saturnTilt, 0.0), SATURN_TILT_MAX);
		}
		break;
	case 'm':
		if (primitiveMode + 1 == sizeof(primitiveTypes)/sizeof(GLenum) && primitiveTypeSphere == false)
			primitiveTypeSphere = true;
		else 
		{
			primitiveTypeSphere = false;
			primitiveMode = (primitiveMode + 1) % (sizeof(primitiveTypes) / sizeof(GLenum));
		}
		break;
	case 'g':
		options[OPT_ANIMATION] = !options[OPT_ANIMATION];
		break;
	case 's':
		options[OPT_RENDER_SATURN] = !options[OPT_RENDER_SATURN];
		break;
	case 'r':
		options[OPT_RENDER_RINGS] = !options[OPT_RENDER_RINGS];
		break;
	case '*':
		saturn_rotate_velocity *= 2;
		break;
	case '/':
		saturn_rotate_velocity /= 2;
		break;
	case 'f':
		options[OPT_CULL_FACE] = !options[OPT_CULL_FACE];
		break;
	case 'l':
		options[OPT_LIGHTING] = !options[OPT_LIGHTING];
		break;
	case 'c':
		options[OPT_CULLING] = !options[OPT_CULLING];
		break;
	case 'n':
		options[OPT_DRAW_NORMALS] = !options[OPT_DRAW_NORMALS];
		break;
	}
	
	updateRenderState();
}

int main(int argc, char **argv)
{
	// Init glut and create window.
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
	glutCreateWindow(WINDOW_TITLE);

	// Set glut defaults.
	glEnable(GL_DEPTH_TEST);

	// Set glut callbacks.
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutReshapeFunc(reshape);
	glutMotionFunc(mouseMove);
	glutPassiveMotionFunc(mouseMove);
	glutMouseFunc(mouseDown);
	glutKeyboardFunc(keyDown);

	// Initialize variables, create scene etc.
	init();

	// Start the main loop (this never returns!).
	glutMainLoop();

	return EXIT_SUCCESS;
}
