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
float MinorDeterminant(int i, int j, mat3 A);


int main(int argc, char* argv[]) {
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    
    LoadTestModel(triangles);
    
    updateCameraAngle(yaw); //initialize camera angle with default yaw

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


void Draw() {

    SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            vec3 dir(x-(SCREEN_WIDTH/2), y-(SCREEN_HEIGHT/2),focalLength);
            Intersection closestIntersection;
            closestIntersection.distance = numeric_limits<float>::max();

            if(ClosestIntersection(cameraPos, dir, triangles, closestIntersection)) {
                //Coloured direct and indirect illumination with shadow
                vec3 colour = triangles[closestIntersection.triangleIndex].color * (DirectLight(closestIntersection)+indirectLight);
                PutPixelSDL( screen, x, y, colour);

                //Coloured direct illumination with shadow
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


bool ClosestIntersection(vec3 start, vec3 dir, const vector<Triangle> &triangles, 
    Intersection &closestIntersection) {

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

        // finished cramer inverse with distance check
        // compute factors needed for determinant
        float a00 = MinorDeterminant(0,0,A);
        float a01 = -MinorDeterminant(1,0,A);
        float a02 = MinorDeterminant(2,0,A);

        // compute determinant
        float det = A[0][0] * a00 + A[1][0] * a01 + A[2][0] * a02;   

        // if determinant is 0, A is not invertible
        if (det == 0 )
            continue;

        // determinant is not 0 => A is invertible => continue computing factors
        float a10 = -MinorDeterminant(0,1,A);
        float a20 = MinorDeterminant(0,2,A);

        // computating the distance t
        vec3 row1 = (1.f/det) * vec3(a00, a10, a20);
        float t = dot(row1, b);

        // t < 0 => no intersection will occur
        if (t < 0)
            continue;

        // compute the rest of the factors
        float a11 = MinorDeterminant(1,1,A);
        float a12 = -MinorDeterminant(2,1,A);
        float a21 = -MinorDeterminant(1,2,A);
        float a22 = MinorDeterminant(2,2,A);

        vec3 col1 = (1.f/det) * vec3(a00, a01, a02);
        vec3 col2 = (1.f/det) * vec3(a10, a11, a12);
        vec3 col3 = (1.f/det) * vec3(a20, a21, a22);
        mat3 invA(col1, col2, col3);

        // compute x = inv(A) * b 
        vec3 x = invA * b;

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


vec3 DirectLight(const Intersection &i) {
    Triangle tri = triangles[i.triangleIndex];
    vec3 e1 = tri.v1-tri.v0;
    vec3 e2 = tri.v2-tri.v0;
    vec3 pos = tri.v0 + i.position.y*e1 + i.position.z*e2;
    vec3 r = lightPos - pos;
    float rsq = glm::dot(r, r);

    float dist = glm::length(r);
    vec3 dir = r * (1.f/dist);

    Intersection j;
    j.distance = numeric_limits<float>::max();
    ClosestIntersection(pos+dir*0.0001f, dir, triangles, j);
    if(j.distance < dist) {
        return vec3(0,0,0);
    }


    vec3 B = lightColour/((float)(4*M_PI*rsq));
    vec3 u_n = tri.normal;
    vec3 u_r = glm::normalize(r);
    vec3 D = B * (max(glm::dot(u_r, u_n), 0.0f));
    return D;
}


float MinorDeterminant(int i, int j, mat3 A) {
    int a1 = i == 0 ? 1 : 0;
    int a2 = i == 2 ? 1 : 2;
    int b1 = j == 0 ? 1 : 0;
    int b2 = j == 2 ? 1 : 2;
    return (A[a1][b1] * A[a2][b2]) - (A[a1][b2] * A[a2][b1]);
}
