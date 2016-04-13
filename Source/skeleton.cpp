#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include <X11/Xlib.h> 
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

struct Light {
    vec3 pos;
    vec3 colour;
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
vec3 cameraPos(0.35f,0.f,-2.f);
float yaw = -M_PI/18.f;
mat3 R;
float cameraSpeed = 0.00002f;
float focus = 0.005f;

//light variables
vector<Light> lights;
float lightSpeed = 0.2f;
vec3 indirectLight = 0.5f * vec3(1,1,1);

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
bool ClosestIntersection(
    vec3 start, 
    vec3 dir, 
    const vector<Triangle> &triangles, 
    Intersection &closestIntersection);
void updateCameraAngle(float angle);
vec3 DirectLight(const Intersection &i, int k);


int main(int argc, char* argv[]) {
    //is necessary for multithreaded access
    XInitThreads();

	screen = InitializeSDL( 2*SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    
    //Load scene triangles
    LoadTestModel(triangles);
    //initialize camera angle with default yaw
    updateCameraAngle(yaw);

    Light ceilingLight;
    ceilingLight.pos = vec3(0, -0.5, -0.7);
    ceilingLight.colour = 14.f * vec3(1,1,1);

    Light fairyLight1;
    fairyLight1.pos = vec3(0.8, -0.6, 0.65);
    fairyLight1.colour = 3.f * vec3(1,0.3,0.7);

    Light fairyLight2;
    fairyLight2.pos = vec3(-0.8, -0.7, 0.65);
    fairyLight2.colour = 3.f * vec3(0.3,1,0);

    Light fairyLight3;
    fairyLight3.pos = vec3(0.8, 0.5, -0.3);
    fairyLight3.colour = 3.f * vec3(1,0.3,0.7);

    Light fairyLight4;
    fairyLight4.pos = vec3(0.8, -0.3, -0.8);
    fairyLight4.colour = 3.f * vec3(1,0.3,0.3);


    lights.push_back(ceilingLight);
    lights.push_back(fairyLight1);
    lights.push_back(fairyLight2);
    lights.push_back(fairyLight3);
    lights.push_back(fairyLight4);

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
    if(keystate[SDLK_UP])
        cameraPos.z += cameraSpeed*dt;
    if(keystate[SDLK_DOWN])
        cameraPos.z -= cameraSpeed*dt; 
    if(keystate[SDLK_LEFT])
        updateCameraAngle(-M_PI/18.f);
    if(keystate[SDLK_RIGHT])
        updateCameraAngle(M_PI/18.f);
    
    //Move light source
    if(keystate[SDLK_w])
        lights[0].pos.z += lightSpeed;
    if(keystate[SDLK_s])
        lights[0].pos.z -= lightSpeed;
    if(keystate[SDLK_d])
        lights[0].pos.x -= lightSpeed;
    if(keystate[SDLK_a])
        lights[0].pos.x += lightSpeed;
    if(keystate[SDLK_q])
        lights[0].pos.y += lightSpeed;
    if(keystate[SDLK_e])
        lights[0].pos.y -= lightSpeed;
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

    float intersections[SCREEN_WIDTH][SCREEN_HEIGHT];

    #pragma omp parallel for
    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            vec3 averageColor(0.0, 0.0, 0.0);
            vec3 averageSectionColor(0.0, 0.0, 0.0); 

            // shoot 9 rays instead of just 1
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
                                *(DirectLight(closestIntersection, 0)
                                + DirectLight(closestIntersection, 1) 
                                + DirectLight(closestIntersection, 2)
                                + DirectLight(closestIntersection, 3)
                                + DirectLight(closestIntersection, 4)    
                                + indirectLight);
                            averageColor += currentColour;
                            intersections[x][y] = closestIntersection.distance;
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
    }

    #pragma omp parallel for
    for(int y = 0; y < SCREEN_HEIGHT; ++y) {
        for(int x = 0; x < SCREEN_WIDTH; ++x) {
            // pixel is in focus, no need to blur out
            if (intersections[x][y] > focus - 0.001f && intersections[x][y] < focus + 0.001f)
                continue;

            vec3 averageColor;

            averageColor =      GetPixelSDL(screen, x-2, y-2) +
                            4.f*GetPixelSDL(screen, x-1, y-2) +
                            7.f*GetPixelSDL(screen, x  , y-2) +
                            4.f*GetPixelSDL(screen, x+1, y-2) +
                                GetPixelSDL(screen, x+2, y-2) +

                            4.f*GetPixelSDL(screen, x-2, y-1) +
                           16.f*GetPixelSDL(screen, x-1, y-1) +
                           26.f*GetPixelSDL(screen, x  , y-1) +
                           16.f*GetPixelSDL(screen, x+1, y-1) +
                            4.f*GetPixelSDL(screen, x+2, y-1) +

                            7.f*GetPixelSDL(screen, x-2, y  ) +
                           26.f*GetPixelSDL(screen, x-1, y  ) +
                           41.f*GetPixelSDL(screen, x  , y  ) +
                           26.f*GetPixelSDL(screen, x+1, y  ) +
                            7.f*GetPixelSDL(screen, x+2, y  ) +

                            4.f*GetPixelSDL(screen, x-2, y+1) +
                           16.f*GetPixelSDL(screen, x-1, y+1) +
                           26.f*GetPixelSDL(screen, x  , y+1) +
                           16.f*GetPixelSDL(screen, x+1, y+1) +
                            4.f*GetPixelSDL(screen, x+2, y+1) +

                                GetPixelSDL(screen, x-2, y+2) +
                            4.f*GetPixelSDL(screen, x-1, y+2) +
                            7.f*GetPixelSDL(screen, x  , y+2) +
                            4.f*GetPixelSDL(screen, x+1, y+2) +
                                GetPixelSDL(screen, x+2, y+2);


            averageColor.x = averageColor.x/273;
            averageColor.y = averageColor.y/273;
            averageColor.z = averageColor.z/273;

            PutPixelSDL( screen, x, y, averageColor);
        }
    }

    ///////////////////////////////////////////////////
    // Roberts operator for edge detection

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

    // edge detection on the grayscale image
    for(int y = 1; y < SCREEN_HEIGHT - 1; ++y) {
        for(int x = 1; x < SCREEN_WIDTH - 1; ++x) {
            float magnitude = abs(GetPixelSDL(screen, x+SCREEN_WIDTH    , y    ).x 
                                - GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y + 1).x)
                            + abs(GetPixelSDL(screen, x+SCREEN_WIDTH    , y + 1).x
                                - GetPixelSDL(screen, x+SCREEN_WIDTH + 1, y    ).x);

            vec3 magVector (magnitude, magnitude, magnitude);
            PutPixelSDL( screen, x+SCREEN_WIDTH, y, magVector);
        }
    }

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
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

        // completed cramer inverse with distance check
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
vec3 DirectLight(const Intersection &i, int k) {
    const Triangle tri = triangles[i.triangleIndex];
    //Triangle coordinate system with origin at v0
    const vec3 e1 = tri.v1-tri.v0;
    const vec3 e2 = tri.v2-tri.v0;
    //Convert triangle-coordinate intersection point to global position
    const vec3 pos = tri.v0 + i.position.y*e1 + i.position.z*e2;
    //Direction from intersection to light source
    const vec3 r = lights[k].pos - pos;
    //distance between intersection and light source
    const float dist = glm::length(r);
    //distance squared between intersection and light source
    const float rsq = glm::dot(r, r);
    //Scaled direction vector with distance to form ray
    //use a small offset (dir*0.0001f) on starting position to prevent intersection with 
    //the originating surface
    const vec3 dir = r/dist ;
    vec3 D;
        
    //Find closest intersection between pixel scene location and the light source
    Intersection j;
    j.distance = numeric_limits<float>::max();

    ClosestIntersection(pos+dir* 0.0001f, dir, triangles, j);

    //If an object is bewteen the surface and the light source, the pixel is in shadow
    if(j.distance < dist) {
        D.x = 0.0;
        D.y = 0.0;
        D.z = 0.0;
        return D;
    }

    // Calculate DirectLight = (lightColour * max(u_r * u_n, 0)) / (4*pi*rsq)
    // Where u_r is the unit direction betwene the surface and the light source,
    // u_n is the unit normals of the triangle, and rsq is the squared distance
    // between the surface and light source.
    // B is the power per area reaching any point in a sphere around the light source,
    const vec3 B = lights[k].colour/((float)(4*M_PI*rsq));
    const vec3 u_n = tri.normal;
    const vec3 u_r = glm::normalize(r);
    D = B * (max(glm::dot(u_r, u_n), 0.0f));  

    return D;
}
