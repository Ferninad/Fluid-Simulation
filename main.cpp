#include "common.h"
#include "cmath"
#include "vector"
#include "OpenSimplexNoise.h"

bool Init();
void CleanUp();
void Run();
void Draw();
void CreateFluid(double ct, double diffusion, double viscosity, int size);
void AddFluid(int x, int y, double dens);
void AddVelocity(int x, int y, double velx, double vely);
void diffuse(int b, vector<vector<double>> &x, vector<vector<double>> &x0);
void project(vector<vector<double>> &velocX, vector<vector<double>> &velocY, vector<vector<double>> &p, vector<vector<double>> &div);
void advect(int b, vector<vector<double>> &d, vector<vector<double>> &d0, vector<vector<double>> &velocX, vector<vector<double>> &velocY);
void lin_solve(int b, vector<vector<double>> &x, vector<vector<double>> &x0, double a, double c);
void set_bnd(int b, vector<vector<double>> &x);
double ScaleNum(double n, double minN, double maxN, double min, double max);
void printDensity();
void printVelocity();
void printS();

SDL_Window *window;
SDL_GLContext glContext;
SDL_Surface *gScreenSurface = nullptr;
SDL_Renderer *renderer = nullptr;
SDL_Rect pos;

int screenWidth = 500;
int screenHeight = 500;
int gridSize = 2;
int N;
int iter = 16;
double dt;
double diff;
double visc;
double t = 0;
double featureSize = 25;
OpenSimplexNoise *noise1 = nullptr;
const double PI = 3.14159265;
int num = 0;

vector<vector<double>> s;
vector<vector<double>> density;
vector<vector<double>> Vx;
vector<vector<double>> Vy;
vector<vector<double>> Vx0;
vector<vector<double>> Vy0;

bool Init()
{
    if (SDL_Init(SDL_INIT_NOPARACHUTE & SDL_INIT_EVERYTHING) != 0)
    {
        SDL_Log("Unable to initialize SDL: %s\n", SDL_GetError());
        return false;
    }
    else
    {
        //Specify OpenGL Version (4.2)
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 4);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);
        SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
        SDL_Log("SDL Initialised");
    }

    //Create Window Instance
    window = SDL_CreateWindow(
        "Game Engine",
        SDL_WINDOWPOS_CENTERED,
        SDL_WINDOWPOS_CENTERED,
        screenWidth,
        screenHeight,   
        SDL_WINDOW_OPENGL);

    //Check that the window was succesfully created
    if (window == NULL)
    {
        //Print error, if null
        printf("Could not create window: %s\n", SDL_GetError());
        return false;
    }
    else{
        gScreenSurface = SDL_GetWindowSurface(window);
        renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
        SDL_Log("Window Successful Generated");
    }
    //Map OpenGL Context to Window
    glContext = SDL_GL_CreateContext(window);

    return true;
}

int main()
{
    //Error Checking/Initialisation
    if (!Init())
    {
        printf("Failed to Initialize");
        return -1;
    }

    // Clear buffer with black background
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);

    //Swap Render Buffers
    SDL_GL_SwapWindow(window);
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);

    Run();

    CleanUp();
    return 0;
}

void CleanUp()
{
    //Free up resources
    SDL_GL_DeleteContext(glContext);
    SDL_DestroyWindow(window);
    SDL_Quit();
}

void Run()
{
    bool gameLoop = true;
    srand(time(NULL));
    long rand1 = rand() * (RAND_MAX + 1) + rand();
    noise1 = new OpenSimplexNoise{rand1};
    
    CreateFluid(0.2, 0, 0.0000001, screenWidth/gridSize); //dt, diff, visc, N
    while (gameLoop)
    {   
        Draw();
        SDL_RenderPresent(renderer);
        pos.x = 0;
        pos.y = 0;
        pos.w = screenWidth;
        pos.h = screenHeight;
        SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
        SDL_RenderFillRect(renderer, &pos);
        
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (event.type == SDL_QUIT)
            {
                gameLoop = false;
            }
            if (event.type == SDL_KEYDOWN)
            {
                switch (event.key.keysym.sym){
                    case SDLK_ESCAPE:
                        gameLoop = false;
                        break;
                    default:
                        break;
                }
            }

            if (event.type == SDL_KEYUP)
            {
                switch (event.key.keysym.sym){
                    default:
                        break;
                }
            }
        }
    }
}

void printDensity(){
    cout << "Density:" << endl;
    for(int i = 0; i < density.size(); i++){
        for(int j = 0; j < density[i].size(); j++){
            cout << density[j][i] << "\t";
        }
        cout << endl;
    }
}
void printS(){
    cout << "S:" << endl;
    for(int i = 0; i < s.size(); i++){
        for(int j = 0; j < s[i].size(); j++){
            cout << s[j][i] << "\t";
        }
        cout << endl;
    }
}

void printVelocity(){
    cout << "Velocity:" << endl;
    for(int i = 0; i < Vx.size(); i++){
        for(int j = 0; j < Vx[i].size(); j++){
            cout << Vx[j][i] << ", " << Vy[j][i] << "\t";
        }
        cout << endl;
    }
}

void Draw(){
    int cx = 0.5 * screenWidth/gridSize;
    int cy = 0.5 * screenHeight/gridSize;
    for(int i = -1; i <= 1; i++){
        for(int j = -1; j <= 1; j++){
            AddFluid(cx+i, cy+j, static_cast<double>(rand() % 101 + 50));
        }
    }
    
    for(int i = 0; i < 2; i++){
        double angle = (*noise1).eval(t/featureSize, 0) * (2*PI) * 2;
        double velx = .2 * cos(angle);
        double vely = .2 * sin(angle);
        t += 0.05;
        AddVelocity(cx, cy, velx, vely);
    }

    diffuse(1, Vx0, Vx);
    diffuse(2, Vy0, Vy);
    
    project(Vx0, Vy0, Vx, Vy);
    advect(1, Vx, Vx0, Vx0, Vy0);
    
    advect(2, Vy, Vy0, Vx0, Vy0);
    
    project(Vx, Vy, Vx0, Vy0);
    diffuse(0, s, density);
    advect(0, density, s, Vx, Vy);

    for(int i = 0; i < density.size(); i++){
        for(int j = 0; j < density[i].size(); j++){
            int color;
            if(density[i][j] > 255)
                color = 255;
            else
                color = ScaleNum(density[i][j], 0, 255, 0, 255);
            pos.x = i * gridSize;
            pos.y = j * gridSize;
            pos.w = gridSize;
            pos.h = gridSize;
            SDL_SetRenderDrawColor(renderer, color, color, color, 255);
            SDL_RenderFillRect(renderer, &pos);
        }
    }
}

void CreateFluid(double ct, double diffusion, double viscosity, int size){
    N = size;
    dt = ct;
    diff = diffusion;
    visc = viscosity;
    
    for(int x = 0; x < N; x++){
        vector<double> temp;
        for(int y = 0; y < N; y++){
            temp.push_back(0);
        }
        s.push_back(temp);
        density.push_back(temp);
        Vx.push_back(temp);
        Vy.push_back(temp);
        Vx0.push_back(temp);
        Vy0.push_back(temp);
    }
}

void AddFluid(int x, int y, double dens){
    density[x][y] += dens;
}

void AddVelocity(int x, int y, double velx, double vely){
    Vx[x][y] += velx;
    Vy[x][y] += vely;
}

void diffuse(int b, vector<vector<double>> &x, vector<vector<double>> &x0){
    double a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}

void lin_solve(int b, vector<vector<double>> &x, vector<vector<double>> &x0, double a, double c){
    double cRecip = 1.0 / c;
    for(int t = 0; t < iter; t++){
        for(int i = 1; i < N - 1; i++){
            for(int j = 1; j < N - 1; j++){
                x[i][j] = (x0[i][j]
                            + a*(x[i+1][j] + x[i-1][j]
                            + x[i][j+1] + x[i][j-1])) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}

void project(vector<vector<double>> &velocX, vector<vector<double>> &velocY, vector<vector<double>> &p, vector<vector<double>> &div){
    for(int i = 1; i < N - 1; i++){
        for(int j = 1; j < N - 1; j++){
            div[i][j] = -0.5*(velocX[i+1][j] - velocX[i-1][j]
            + velocY[i][j+1] - velocY[i][j-1])/N;
            p[i][j] = 0;
        }
    }

    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);

    for(int i = 1; i < N - 1; i++){
        for(int j = 1; j < N - 1; j++){
            velocX[i][j] -= 0.5 * (p[i+1][j] - p[i-1][j]) * N;
            velocY[i][j] -= 0.5 * (p[i][j+1] - p[i][j-1]) * N;
        }
    }
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

void advect(int b, vector<vector<double>> &d, vector<vector<double>> &d0, vector<vector<double>> &velocX, vector<vector<double>> &velocY){
    double i0, i1, j0, j1;
  
    double dtx = dt * (N - 2);
    double dty = dt * (N - 2);
    
    double s0, s1, t0, t1;
    double tmp1, tmp2, tmp3, x, y;
  
    double Nfloat = N;
    double ifloat, jfloat;
    int i, j, k;
    

    for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++){
        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++){
            tmp1 = dtx * velocX[i][j];
            tmp2 = dty * velocY[i][j];
    
            x = ifloat - tmp1;
            y = jfloat - tmp2;
    
            if(x < 0.5) x = 0.5;
            if(x > Nfloat + 0.5) x = Nfloat + 0.5;
            i0 = floor(x);
            i1 = i0 + 1.0;
            if(y < 0.5) y = 0.5; 
            if(y > Nfloat + 0.5) y = Nfloat + 0.5;
            j0 = floor(y);
            j1 = j0 + 1.0;
            
            s1 = x - i0; 
            s0 = 1.0 - s1; 
            t1 = y - j0; 
            t0 = 1.0 - t1;
        
            int i0i = i0;
            int i1i = i1;
            int j0i = j0;
            int j1i = j1;
            
            d[i][j] = 
			s0 * (t0 * d0[i0i][j0i] + t1 * d0[i0i][j1i]) +
			s1 * (t0 * d0[i1i][j0i] + t1 * d0[i1i][j1i]);
        }
    }
    
    set_bnd(b, d);
}

void set_bnd(int b, vector<vector<double>> &x){
    for(int i = 1; i < N - 1; i++){
		x[i][0] = b == 2 ? -x[i][1] : x[i][1];
		x[i][N-1] = b == 2 ? -x[i][N-2] : x[i][N-2];
		}
	for(int j = 1; j < N - 1; j++){
		x[0][j] = b == 1 ? -x[1][j] : x[1][j];
		x[N-1][j] = b == 1 ? -x[N-2][j] : x[N-2][j];
	}

	x[0][0] = 0.5 * (x[1][0] + x[0][1]);
	x[0][N-1] = 0.5 * (x[1][N-1] + x[0][N-2]);
	x[N-1][0] = 0.5 * (x[N-2][0] + x[N-1][1]);
    x[N-1][N-1] = 0.5 * (x[N-2][N-1] + x[N-1][N-2]);
}

double ScaleNum(double n, double minN, double maxN, double min, double max){
    return (((n - minN) / (maxN - minN)) * (max - min)) + min;
}