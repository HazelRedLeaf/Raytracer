#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

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
int t;
vector<Triangle> triangles;

//camera variables
float focalLength = 250.f;
vec3 cameraPos(0.f,0.f,-1.5f);
float yaw = -M_PI/18.f;
mat3 R;
float cameraSpeed = 0.2f;

//light variables
vec3 lightPos(0, -0.5, -0.7);
vec3 lightColour = 14.f * vec3(1,1,1);
float lightSpeed = 0.2f;
vec3 indirectLight = 0.5f*vec3(1,1,1);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, Intersection &closestIntersection);
void updateCameraAngle(float angle);
vec3 DirectLight(const Intersection &i);
//float MinorDeterminant(int i, int j, mat3 A);

int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    
    LoadTestModel(triangles);
    //initialize camera angle with default yaw
    updateCameraAngle(yaw);

	while( NoQuitMessageSDL() )
	{
		Update();
		Draw();
	}

	SDL_SaveBMP( screen, "screenshot.bmp" );

	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2-t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

    Uint8* keystate = SDL_GetKeyState(0);
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
    R = mat3(vec3(cos(angle),0,sin(angle)),vec3(0,1,0),vec3(-sin(angle),0,cos(angle)));
    cameraPos = R * cameraPos;
}

void Draw()
{

    SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            vec3 dir(x-(SCREEN_WIDTH/2), y-(SCREEN_HEIGHT/2),focalLength);
            Intersection closestIntersection;
            closestIntersection.distance = numeric_limits<float>::max();
            if(ClosestIntersection(cameraPos, dir, triangles, closestIntersection)) {
                //Coloured direct and indirect illumination
                vec3 colour = triangles[closestIntersection.triangleIndex].color * (DirectLight(closestIntersection)+indirectLight);
                PutPixelSDL( screen, x, y, colour);
                //Coloured direct illumination
                //vec3 colour = triangles[closestIntersection.triangleIndex].color * DirectLight(closestIntersection);
                //PutPixelSDL( screen, x, y, colour);
                //Greyscale direct illumination
                //PutPixelSDL( screen, x, y, DirectLight(closestIntersection));
                //Colour, no light
                //PutPixelSDL( screen, x, y, triangles[closestIntersection.triangleIndex].color);
            }
            else {
                PutPixelSDL( screen, x, y, vec3(0,0,0));
            }
        }
    }

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}
/*
float MinorDeterminant(int i, int j, mat3 A) {
    int a1 = i == 0 ? 1 : 0;
    int a2 = i == 2 ? 1 : 2;
    int b1 = j == 0 ? 1 : 0;
    int b2 = j == 2 ? 1 : 2;
    return (A[a1][b1] * A[a2][b2]) - (A[a1][b2] * A[a2][b1]);
}*/

bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, Intersection &closestIntersection) {
    bool foundIntersect = false;

    for(unsigned int i=0; i < triangles.size(); i++) {
        Triangle triangle = triangles[i];
        vec3 v0 = triangle.v0;
        vec3 v1 = triangle.v1;
        vec3 v2 = triangle.v2;
        
        vec3 e1 = v1-v0;
        vec3 e2 = v2-v0;
        vec3 b = start-v0;

        mat3 A(-dir,e1,e2);
        vec3 x = glm::inverse(A)*b;
        //unfinished cramer inverse with distance check
        /*float a11 = MinorDeterminant(0,0,A);
        float a12 = -MinorDeterminant(1,0,A);
        float a13 = MinorDeterminant(2,0,A);
        
        float det = A[0][0] * a11 - A[1][0] * a12 + A[2][0] * a13;        

        float a21 = -MinorDeterminant(0,1,A);
        float a31 = MinorDeterminant(0,2,A);
         
        vec3 row1 = vec3(a11,a21,a31) * (1.f/det);

        float t = b.x*row1.x + b.y * row1.y + b.z * row1.z;
        if(t < 0) {
            continue;  
        }
        
        float a22 = MinorDeterminant(1,1,A);
        float a23 = -MinorDeterminant(2,1,A);
        float a32 = -MinorDeterminant(1,2,A);
        float a33 = MinorDeterminant(2,2,A);
        vec3 row2(a12,a22,a32);
        vec3 row3(a13,a23,a33);
        
        vec3 x = (mat3(row1,row2,row3) * (1.f/det)) * b;*/

        if(0 < x.y && 0 < x.z && (x.y + x.z) < 1 && 0 <= x.x) {
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

vec3 DirectLight(const Intersection &i) {
    Triangle tri = triangles[i.triangleIndex];
    vec3 e1 = tri.v1-tri.v0;
    vec3 e2 = tri.v2-tri.v0;
    vec3 r = tri.v0 + i.position.y*e1 + i.position.z*e2;
    r = lightPos - r;
    float rsq = glm::dot(r, r);
    vec3 B = lightColour/((float)(4*M_PI*rsq));
    vec3 u_n = tri.normal;
    vec3 u_r = glm::normalize(r);
    vec3 D = B * (max(glm::dot(u_r, u_n), 0.0f));
    return D;
}
