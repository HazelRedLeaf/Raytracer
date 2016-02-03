#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModel.h"

using namespace std;
using glm::vec3;
using glm::mat3;

/* ----------------------------------------------------------------------------*/
/* GLOBAL VARIABLES                                                            */

const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<vec3> stars(1000);
const float velocity = 0.0001f;
/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

void Update();
void Draw();
void Interpolate(float a, float b, vector<float> &result);
void Interpolate(vec3 a, vec3 b, vector<vec3> &result);


int main( int argc, char* argv[] )
{
	screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT );
	t = SDL_GetTicks();	// Set start value for timer.
    
    
    for(int i=0; i < stars.size(); i++) {
        float x = -1 + 2*float(rand()) / float(RAND_MAX);
        float y = -1 + 2*float(rand()) / float(RAND_MAX);
        float z = float(rand()) / float(RAND_MAX);
        stars[i] = vec3(x,y,z);
    }

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
    for(int i=0; i < stars.size(); i++) {
        stars[i].z -= velocity * dt;
        if(stars[i].z <= 0)
            stars[i].z += 1;
        if(stars[i].z > 1)
            stars[i].z -= 1;
    }
}

void Draw()
{

    SDL_FillRect(screen, 0, 0);

	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

    float f = SCREEN_HEIGHT/2;
    float g = SCREEN_WIDTH/2;

    for(int i=0; i<stars.size();i++) {
        float u =  f * stars[i].x/stars[i].z + g;
        float v = f * stars[i].y/stars[i].z + f;
        vec3 colour = 0.2f * vec3(1,1,1)/(stars[i].z*stars[i].z);
        PutPixelSDL( screen, (int) u, (int) v, colour);
    }

	if( SDL_MUSTLOCK(screen) )
		SDL_UnlockSurface(screen);

	SDL_UpdateRect( screen, 0, 0, 0, 0 );
}

void Interpolate(float a, float b, vector<float> &result) {
	// return average of a and b if result of size 1
	if(result.size() == 1) {
		result[0] = a+(b-a)/2;
	}
	// do nothing if result of size 0
	else if(result.size() == 0) {
		return;
	}
	float step_size = (b-a)/(result.size()-1);
	//cout << "Hello";
	for(int i=0; i<result.size(); i++) {
		result[i] = a+i*step_size;
	}
}

void Interpolate(vec3 a, vec3 b, vector<vec3> &result) {
	// return average of a and b if result of size 1
    vec3 d = b - a;
    
	if(result.size() == 1) {
		result[0] = vec3(a.x+d.x/2,a.y+d.y/2,a.z+d.z/2);
	}
	// do nothing if result of size 0
	else if(result.size() == 0) {
		return;
	}
	vec3 step_size = d/(float)(result.size()-1);

	for(int i=0; i<result.size(); i++) {
		result[i] = a+(float)i*step_size;
	}
}
