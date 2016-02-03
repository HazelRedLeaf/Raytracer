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

    vector<vec3> stars(5);
    
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
}

void Draw()
{
	if( SDL_MUSTLOCK(screen) )
		SDL_LockSurface(screen);

	vec3 top_left(1,0,0);
    vec3 top_right(0,0,1);
    vec3 bottom_right(0,1,0);
    vec3 bottom_left(1,1,0);

    vector<vec3> leftSide(SCREEN_HEIGHT);
    vector<vec3> rightSide(SCREEN_HEIGHT);
    Interpolate(top_left, bottom_left, leftSide);
    Interpolate(top_right, bottom_right, rightSide);

	for( int y=0; y<SCREEN_HEIGHT; ++y )
	{
        vector<vec3> row(SCREEN_WIDTH);
        Interpolate(leftSide[y], rightSide[y], row);
		for( int x=0; x<SCREEN_WIDTH; ++x )
		{
			PutPixelSDL( screen, x, y, row[x] );
		}
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
