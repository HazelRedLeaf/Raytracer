#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <X11/Xlib.h> 
#include "SDLauxiliary.h"
#include "TestModel.h"
// header inclusion for linux lab machines
#include "/usr/include/opencv2/objdetect/objdetect.hpp"
#include "/usr/include/opencv2/opencv.hpp"
#include "/usr/include/opencv2/core/core.hpp"
#include "/usr/include/opencv2/highgui/highgui.hpp"
#include "/usr/include/opencv2/imgproc/imgproc.hpp"

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* STRUCTS                                                                     */
struct Intersection {
    vec3 position;
    float distance;
    int triangleIndex;
};

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
SDL_Surface* edgeScreen;
int t;
vector<Triangle> triangles;

//camera variables
float focalLength = 300.f;
vec3 cameraPos(0.2f,0.f,-2.f);
float yaw = -M_PI/18.f;
mat3 R;
float cameraSpeed = 0.2f;

//light variables
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColour = 14.f * vec3(1,1,1);
float lightSpeed = 0.2f;
vec3 indirectLight = 0.5f * vec3(1,1,1);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, 
        Intersection &closestIntersection);
void updateCameraAngle(float angle);
vec3 DirectLight(const Intersection &i);


int main(int argc, char* argv[]) {
    //is necessary for multithreaded access
    XInitThreads();

	screen = InitializeSDL( 2*SCREEN_WIDTH, SCREEN_HEIGHT );
    //edgeScreen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    
    //Load scene triangles
    LoadTestModel(triangles);
    //initialize camera angle with default yaw
    updateCameraAngle(yaw); 

	while(NoQuitMessageSDL()) {
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	return 0;
}


void Update() {
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

    Uint8* keystate = SDL_GetKeyState(0);
    //Move camera
    if(keystate[SDLK_UP]) {
        cameraPos.z += cameraSpeed;
    } 
    if(keystate[SDLK_DOWN]) {
        cameraPos.z -= cameraSpeed;
    } 
    if(keystate[SDLK_LEFT]) {
        cameraPos.x -= cameraSpeed;
        updateCameraAngle(-M_PI/18.f);
    } 
    if(keystate[SDLK_RIGHT]) {
        cameraPos.x += cameraSpeed;
        updateCameraAngle(M_PI/18.f);
    } 
    //Move light source
    if(keystate[SDLK_w]) {
        lightPos.z += lightSpeed;
    }
    if(keystate[SDLK_s]){
        lightPos.z -= lightSpeed;
    }
    if(keystate[SDLK_d]){
        lightPos.x -= lightSpeed;
    }
    if(keystate[SDLK_a]){
        lightPos.x += lightSpeed;
    }
    if(keystate[SDLK_q]){
        lightPos.y += lightSpeed;
    }
    if(keystate[SDLK_e]){
        lightPos.y -= lightSpeed;
    }
}


void updateCameraAngle(float angle) {
    yaw += angle;
    //update rotation matrix with angle
    R = mat3(vec3( cos(angle), 0, sin(angle)),
             vec3(      0,     1,   0    ),
             vec3(-sin(angle), 0, cos(angle)));
    //update camera position with rotation matrix
    cameraPos = R * cameraPos;
}


void Draw() {

    SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    #pragma omp parallel for
    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            vec3 averageColor(0.0,0.0,0.0);

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    //ray direction from current pixel
                    vec3 dir(x-(SCREEN_WIDTH/2)+i, y-(SCREEN_HEIGHT/2)+j,focalLength);
                    //Find the closest intersected triangle from the current pixel/camera position
                    Intersection closestIntersection;
                    closestIntersection.distance = numeric_limits<float>::max();
                    if(ClosestIntersection(cameraPos, dir, triangles, closestIntersection)) {
                        //use the intersected triangle to find the pixel's colour and illumination/shadow.
                        //Coloured direct and indirect illumination with shadow
                        vec3 currentColour = triangles[closestIntersection.triangleIndex].color 
                            * (DirectLight(closestIntersection)+indirectLight);
                        averageColor += currentColour;
                        //PutPixelSDL( screen, x, y, colour);
                    }
                    // No intersection found (eg outside of scene bounds) so colour pixel black
                    else {
                        PutPixelSDL( screen, x, y, vec3(0,0,0));
                    }
                }
            }

            averageColor.x = averageColor.x/9;
            averageColor.y = averageColor.y/9;
            averageColor.z = averageColor.z/9;

            PutPixelSDL( screen, x, y, averageColor);
                    
        }
        //SDL_UpdateRect( screen, 0, 0, 0, 0 );
    }

    ///////////////////////////////////////////////////
    // Sobel edge detection attempt

    // kernels
    mat3 dXkernel = mat3(vec3( -1, -2, -1),
                         vec3(  0,  0,  0),
                         vec3(  1,  2,  1)); 

    mat3 dYkernel = mat3(vec3( -1,  0,  1),
                         vec3( -2,  0,  2),
                         vec3( -1,  0,  1));

    // Grayscale representation of current frame
    #pragma omp parallel for
    for(int y = 1; y < SCREEN_HEIGHT - 1; ++y) {
        for(int x = 1; x < SCREEN_WIDTH - 1; ++x) {
            vec3 color = GetPixelSDL(screen, x, y);
            float graycolor = (color.x + color.y + color.z) / 3;
            vec3 grayColorVec (graycolor, graycolor, graycolor);
            PutPixelSDL( screen, x+SCREEN_WIDTH, y, grayColorVec);
        }
    }

    // convolution
    //#pragma omp parallel for
    for(int y = 1; y < SCREEN_HEIGHT - 1; ++y) {
        for(int x = 1; x < SCREEN_WIDTH - 1; ++x) {
            float magnitude = abs(GetPixelSDL(screen, x+SCREEN_WIDTH    , y    ).x 
                                - GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x)
                            + abs(GetPixelSDL(screen, x+SCREEN_WIDTH    , y + 1).x
                                - GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y    ).x);

            // float pixelX = (dXkernel[0][0] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x +
            //                 dXkernel[0][1] * GetPixelSDL(screen, x+SCREEN_WIDTH    , y + 1).x +
            //                 dXkernel[0][2] * GetPixelSDL(screen, x+SCREEN_WIDTH - 1, y + 1).x +
            //                 dXkernel[1][0] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y    ).x +
            //                 dXkernel[1][1] * GetPixelSDL(screen, x+SCREEN_WIDTH,     y    ).x +
            //                 dXkernel[1][2] * GetPixelSDL(screen, x+SCREEN_WIDTH - 1, y    ).x +
            //                 dXkernel[2][0] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x +
            //                 dXkernel[2][1] * GetPixelSDL(screen, x+SCREEN_WIDTH,     y - 1).x +
            //                 dXkernel[2][2] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x );

            // float pixelY = (dYkernel[0][0] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x +
            //                 dYkernel[0][1] * GetPixelSDL(screen, x+SCREEN_WIDTH    , y + 1).x +
            //                 dYkernel[0][2] * GetPixelSDL(screen, x+SCREEN_WIDTH - 1, y + 1).x +
            //                 dYkernel[1][0] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y    ).x +
            //                 dYkernel[1][1] * GetPixelSDL(screen, x+SCREEN_WIDTH,     y    ).x +
            //                 dYkernel[1][2] * GetPixelSDL(screen, x+SCREEN_WIDTH - 1, y    ).x +
            //                 dYkernel[2][0] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x +
            //                 dYkernel[2][1] * GetPixelSDL(screen, x+SCREEN_WIDTH,     y - 1).x +
            //                 dYkernel[2][2] * GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x );

            // float magnitude = sqrt(pixelX * pixelX + pixelY * pixelY);
            // magnitude = magnitude * 255.0 / (4 * 255 * sqrt(2));
            vec3 magVector (magnitude, magnitude, magnitude);
            PutPixelSDL( screen, x+SCREEN_WIDTH, y, magVector);
        }
    }

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
    //SDL_UpdateRect( edgeScreen, 0, 0, 0, 0 );
}


bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, 
    Intersection &closestIntersection) {

    bool foundIntersect = false;

    //Iterate through all triangles in scene to find where the ray intersects each (if at all)
    //and to find the closest intersection to the start position
    for(unsigned int i=0; i < triangles.size(); i++) {
        const Triangle triangle = triangles[i];
        const vec3 v0 = triangle.v0;
        const vec3 v1 = triangle.v1;
        const vec3 v2 = triangle.v2;
        
        //create triangle coordinate system (with origin at v0)
        const vec3 e1 = v1-v0;
        const vec3 e2 = v2-v0;
        const vec3 b = start-v0;

        //compose matrix A of the negative direction vector and the triangle edges
        const mat3 A(-dir,e1,e2);
        // finished cramer inverse with distance check
        // compute factors needed for determinant
        float a00 = (A[1][1] * A[2][2]) - (A[1][2] * A[2][1]);
        float a01 = -((A[0][1] * A[2][2]) - (A[0][2] * A[2][1]));
        float a02 = (A[0][1] * A[1][2]) - (A[0][2] * A[1][1]);

        // compute determinant
        const float det = A[0][0] * a00 + A[1][0] * a01 + A[2][0] * a02;   
        const float invDet = 1.f/det;

        // if determinant is 0, A is not invertible
        if (det == 0 )
            continue;

        a00 *= invDet;
        a01 *= invDet;
        a02 *= invDet;
        // determinant is not 0 => A is invertible => continue computing factors
        const float a10 = -((A[1][0] * A[2][2]) - (A[1][2] * A[2][0])) * invDet;
        const float a20 = ((A[1][0] * A[2][1]) - (A[1][1] * A[2][0])) * invDet;

        // computating the distance t
        const vec3 row1 = vec3(a00, a10, a20);
        const float t = dot(row1, b);

        // t < 0 => no intersection will occur
        if (t < 0)
            continue;

        // compute the rest of the factors
        const float a11 = ((A[0][0] * A[2][2]) - (A[0][2] * A[2][0])) * invDet;
        const float a12 = -((A[0][0] * A[1][2]) - (A[0][2] * A[1][0])) * invDet;
        const float a21 = -((A[0][0] * A[2][1]) - (A[0][1] * A[2][0])) * invDet;
        const float a22 = ((A[0][0] * A[1][1]) - (A[0][1] * A[1][0])) * invDet;

        //Put inverse results together
        const vec3 col1 = vec3(a00, a01, a02);
        const vec3 col2 = vec3(a10, a11, a12);
        const vec3 col3 = vec3(a20, a21, a22);
        const mat3 invA = mat3(col1, col2, col3);

        // compute x = inv(A) * b 
        const vec3 x = invA * b;
        //vec3 x = glm::inverse(A)*b;
        if(0 < x.y && 0 < x.z && (x.y + x.z) < 1 && 0 < x.x) {
            if(closestIntersection.distance > x.x) {
                closestIntersection.triangleIndex = i;
                closestIntersection.position = x;
                closestIntersection.distance = x.x;
                foundIntersect = true;
            }
        }
    }
    return foundIntersect;
}

//Calculates the direct light for a pixel given its closest intersection point
vec3 DirectLight(const Intersection &i) {
    const Triangle tri = triangles[i.triangleIndex];
    //Triangle coordinate system with origin at v0
    const vec3 e1 = tri.v1-tri.v0;
    const vec3 e2 = tri.v2-tri.v0;
    //Convert triangle-coordinate intersection point to global position
    const vec3 pos = tri.v0 + i.position.y*e1 + i.position.z*e2;
    //Direction from intersection to light source
    const vec3 r = lightPos - pos;
    //distance between intersection and light source
    const float dist = glm::length(r);
    //distance squared between intersection and light source
    const float rsq = glm::dot(r, r);
    //Scaled direction vector with distance to form ray
    const vec3 dir = r * (1.f/dist);

    //Find closest intersection between pixel scene location and the light source
    Intersection j;
    j.distance = numeric_limits<float>::max();
    //use a small offset (dir*0.0001f) on starting position to prevent intersection with 
    //the originating surface
    ClosestIntersection(pos+dir*0.0001f, dir, triangles, j);
    //If an object is bewteen the surface and the light source, the pixel is in shadow
    if(j.distance < dist) {
        return vec3(0,0,0);
    }

    //Calculate DirectLight = (lightColour * max(u_r * u_n, 0)) / (4*pi*rsq)
    //Where u_r is the unit direction betwene the surface and the light source,
    //u_n is the unit normals of the triangle, and rsq is the squared distance
    //between the surface and light source.
    //B is the power per area reaching any point in a sphere around the light source,
    const vec3 B = lightColour/((float)(4*M_PI*rsq));
    const vec3 u_n = tri.normal;
    const vec3 u_r = glm::normalize(r);
    const vec3 D = B * (max(glm::dot(u_r, u_n), 0.0f));
    return D;
}
