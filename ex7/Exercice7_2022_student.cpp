#include "ConfigFile.tpp"
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>


#define M_PI 3.14159265358979323846

enum class BoundaryType
{
    Forcing,
    Dirichlet,
    Exit,
    Neumann,
};

BoundaryType
boundaryConditionFromString(const std::string& s)
{
    if (s == "Forcing")
        return BoundaryType::Forcing;
    else if (s == "Dirichlet")
        return BoundaryType::Dirichlet;
    else if (s == "Exit")
        return BoundaryType::Exit;
    else if (s == "Neumann")
        return BoundaryType::Neumann;
    else
        throw std::runtime_error("unrecognised BoundaryType string");
}

struct BoundaryConditions
{
    BoundaryType left, right, up, down;
};

enum class InitialCondition
{
    Zero,
    FinitStatic,
    FinitRight,
    FinitLeft,
    Eigenmode
};

InitialCondition
initialConditionFromString(const std::string& s)
{
    if (s == "Zero")
        return InitialCondition::Zero;
    else if (s == "FinitStatic")
        return InitialCondition::FinitStatic;
    else if (s == "FinitRight")
        return InitialCondition::FinitRight;
    else if (s == "FinitLeft")
        return InitialCondition::FinitLeft;
    else if (s == "Eigenmode")
        return InitialCondition::Eigenmode;
    else
        throw std::runtime_error("unrecognised InitialCondition string");
}

using namespace std;

void
writeData(const string& filepath,
          const double time,
          const vector<double>& x,
          const vector<double>& y,
          const vector<vector<double>>& heightField)
{
    const int xDim = x.size();
    const int yDim = y.size();

    if (xDim != heightField.size())
        throw std::runtime_error("size mismatch: x and heightfield different dimensions");

    ofstream os(filepath);
    os.precision(15);
    // First store current simulation time
    os << time << "\n";

    // Then store all data
    for (unsigned int i = 0; i < xDim; ++i) {
        if (yDim != heightField[i].size())
            throw std::runtime_error("size mismatch: y and heightfield different dimensions");
        for (int j = 0; j < yDim; ++j) {
            const double f = heightField[i][j];
            os << x[i] << " " << y[j] << " " << f << "\n";
        }
    }
}

// Replaces '#' in the filepath with the number of the current frame
string
getFilepathForFrame(int frame, const string& filepath)
{
    const int hashPos = filepath.find('#');
    if (hashPos == string::npos) return filepath;
    std::stringstream ss;
    ss << filepath.substr(0, hashPos);
    ss << to_string(frame);
    ss << filepath.substr(hashPos + 1);
    return ss.str();
}

// The 'std::function<double(double, double)>' is a way to describe an argument that is a function
// with the following signature
//    double myFunction(double argument1, double argument2);
// and in our case makes it possible to use the u2 function to be given as an argument: this reduces
// the amount of parameters we will need to pass around and encapsulates this function better.
using WaveSpeedFunction = std::function<double(double, double)>;

WaveSpeedFunction
u2FunctionGenerator(const std::string& propagationFunctionName,
                    const double u0,
                    const double g,
                    const double h0,
                    const double h1,
                    const double a,
                    const double b,
                    const double Ly)
{
    // Note that the output type here is a function
    if (propagationFunctionName == "Const") {
        return [=](const double x, const double y) -> double {
            // TODO: Implement a constant u2 function
            return u0*u0;
        };
    } else if (propagationFunctionName == "Case1") {
        return [=](const double x, const double y) -> double {
            // TODO: Implement the function defined by "case 1" in the problem sheet
            double h(h0);
            if(a<x && x<b){
                h-=(h0-h1)/2.0*(1.0-cos(M_PI*(x-a)/(b-a)));
            }else if(b<=x){
                h=h1;
            }
            
            return g*h;
        };
    } else if (propagationFunctionName == "Case2") {
        return [=](const double x, const double y) -> double {
            // TODO: Implement the function defined by "case 2" in the problem sheet
            double h(h0);
            if(a<x&&x<b){
                h-=(h0-h1)/2.0*(1.0-cos(M_PI*(x-a)/(b-a)))*sqrt(sin(M_PI*y/Ly));
            }else if(b<=x){
                h-=(h0-h1)*sqrt(sin(M_PI*y/Ly));
            }
            return g*h;
        };
    } else
        throw std::runtime_error("unrecognised name for PropagationFunction");
}

double
calculateDt(const double betaCFL,
            const vector<double>& x,
            const vector<double>& y,
            const WaveSpeedFunction u2)
{
    if (x.empty() || y.empty()) throw std::runtime_error("vector must be larger than 0");
    double u2max(u2(x[0],y[0]));
    for(auto& xi:x){
        for(auto& yj:y){
            if(u2(xi,yj)>u2max){
                u2max=u2(xi,yj);
            }
        }
    }
    // TODO: Implement time-step CFL condition
    return betaCFL/(sqrt(u2max))/sqrt(1.0/pow(x[1]-x[0],2)+1.0/pow(y[1]-y[0],2));
}


void
setInitialConditions(vector<vector<double>>& height,
                     vector<vector<double>>& heightPrev,
                     InitialCondition initialCondition,
                     const vector<double>& x,
                     const vector<double>& y,
                     const double dt,
                     const WaveSpeedFunction& u2,
                     const double f0,
                     const double xa,
                     const double xb,
                     const double ya,
                     const double yb,
                     const double Lx,
                     const double Ly,
                     const double eigenmodeAmplitude,
                     const double eigenmodeM,
                     const double eigenmodeN)
{
    // TODO: Implement initial conditions
    switch (initialCondition) {
        case InitialCondition::Zero:
            for(size_t i(0);i<height.size();++i){
                for(size_t j(0);j<height[i].size();++j){
                    height[i][j]=0;
                    heightPrev[i][j]=0;
                }
            }
            break;
        // Note: Finit refers to the function specified for ex. 7.1f.
        case InitialCondition::FinitStatic:
            for(size_t i(0);i<height.size();++i){
                for(size_t j(0);j<height[i].size();++j){
                    if(x[i]<=xa || xb<=x[i] || y[j]<=ya || yb<=y[j]){
                        height[i][j]=0;
                        heightPrev[i][j]=0;
                    }else{
                        height[i][j]=f0*(1-cos(2*M_PI*(x[i]-xa)/(xb-xa)))*(1-cos(2*M_PI*(y[j]-ya)/(yb-ya)));
                    }
                    heightPrev[i][j]=height[i][j];
                }
            }
            break;
        // Note: Finit refers to the function specified for ex. 7.1f.
             case InitialCondition::FinitRight:
            for(size_t i(0);i<height.size();++i){
                for(size_t j(0);j<height[i].size();++j){
                    if(x[i]<=xa || xb<=x[i] || y[j]<=ya || yb<=y[j]){
                        height[i][j]=0;
                        heightPrev[i][j]=0;
                    }else{
                        height[i][j]=f0*(1-cos(2*M_PI*(x[i]-xa)/(xb-xa)))*(1-cos(2*M_PI*(y[j]-ya)/(yb-ya)));
                    }
                    unsigned int l(round(sqrt(u2(x[i],y[j])*dt*height.size()/(1.0*Lx))));
                    if(i+l<height.size()){
                        heightPrev[i][j]=height[i+l][j];
                    }
                }
            }
            break;
        case InitialCondition::FinitLeft:
            for(size_t i(0);i<height.size();++i){
                for(size_t j(0);j<height[i].size();++j){
                    if(x[i]<=xa || xb<=x[i] || y[j]<=ya || yb<=y[j]){
                        height[i][j]=0;
                        heightPrev[i][j]=0;
                    }else{
                        height[i][j]=f0*(1-cos(2*M_PI*(x[i]-xa)/(xb-xa)))*(1-cos(2*M_PI*(y[j]-ya)/(yb-ya)));
                    }
                    unsigned int l(round(sqrt(u2(x[i],y[j])*dt*height.size()/(1.0*Lx))));
                    if(i-l>0){
                        heightPrev[i][j]=height[i-l][j];
                    }
                }
            }
            break;
        case InitialCondition::Eigenmode:
			for(size_t i(0); i<height.size();++i){
                for(size_t j(0); j<height[i].size();++i){
                    double km(eigenmodeM*M_PI/Lx);
                    double kn(eigenmodeN*M_PI/Ly);
                    height[i][j]=eigenmodeAmplitude*cos(km*x[i]+kn*y[j]);
                    heightPrev[i][j]=height[i][j]*cos(-dt*sqrt(pow(eigenmodeM*M_PI/Lx,2.0)+pow(eigenmodeN*M_PI/Ly,2.0))*sqrt(u2(x[i],y[j])));
                }
            }
            break;
        default:
            height = vector<vector<double>>(x.size(), vector<double>(y.size()));
            heightPrev = height;
    }
}

void
setBoundaryConditions(const BoundaryConditions boundaryConditions,
                      const double time,
                      const double dt,
                      vector<vector<double>>& heightNext,
                      const vector<vector<double>>& height,
                      const vector<double>& x,
                      const vector<double>& y,
                      const double A,
                      const double omega,
                      const WaveSpeedFunction u2)
{
	///fonction probablement améliorable.
    // TODO: Implement boundary conditions.
    
    ///bord gauche/left
    if(boundaryConditions.left==BoundaryType::Neumann){
		heightNext[0]=height[0];
	}
	if(boundaryConditions.left==BoundaryType::Dirichlet){
		for(auto& f:heightNext[0]){
			f=0.0;
		}
	}
	if(boundaryConditions.left==BoundaryType::Forcing){
		for(auto& f : heightNext[0]){
				f=A*sin(omega*(time+dt));
				//cout<<"DEBUG:: Forcing condition actived : f = "<<f<<endl;
		}
		//cout<<"A"<<endl;
	}
	
	///bord droite/right
    if(boundaryConditions.right==BoundaryType::Neumann){
		heightNext.back()=heightNext[heightNext.size()-2];
		//cout<<"B"<<endl;
		
	}
	if(boundaryConditions.right==BoundaryType::Dirichlet){
		for(auto& f : heightNext.back())
		{
			f=0.0;
		}		
	}
	if(boundaryConditions.right==BoundaryType::Exit){
		unsigned int x_ = heightNext[0].size()-1;
		for(unsigned int i(0);i<heightNext.back().size()-1;i++){
			double hx = x[i+1]-x[i];
			//heightNext.back()[i]=height.back()[i]	+sqrt(u2(x_,i))*dt/hx*(height[x_][i]-height[x_+1][i]);
			heightNext.back()[i]=height.back()[i]	+sqrt(u2(x_,i))*dt/hx*(height[x_-1][i]-height.back()[i]);
		}
	}
    
    ///bord bas
    if(boundaryConditions.down==BoundaryType::Neumann){
		for(unsigned int i(0); i<heightNext.size();i++){
			heightNext[i][0]=heightNext[i][1];
		}
		//cout<<"C"<<endl;
	}
	if(boundaryConditions.down==BoundaryType::Dirichlet){
		for(unsigned int i(0); i<heightNext[0].size();i++){
			//heightNext[0][i]=0.0;
		}		
		
	}
    ///bord haut
    if(boundaryConditions.up==BoundaryType::Neumann){
		for(unsigned int i(0); i<heightNext.size();i++){
			heightNext[i][heightNext.back().size()-1]=heightNext[i][heightNext.back().size()-2];
		}
		//cout<<"D"<<endl;
		
	}
	if(boundaryConditions.up==BoundaryType::Dirichlet){
		for(unsigned int i(0); i<heightNext[0].size();i++){
			heightNext.back()[i]=0.0;
		}
	}
}

void
advanceHeightField(const double time,
                   const double dt,
                   vector<vector<double>>& heightNext,
                   const vector<vector<double>>& height,
                   const vector<vector<double>>& heightPrev,
                   const vector<double>& x,
                   const vector<double>& y,
                   const WaveSpeedFunction u2)
{
    // TODO: Implement explicit time integration function.
    for (unsigned int i = 1; i < x.size() - 1; ++i) {
        for (unsigned int j = 1; j < y.size() - 1; ++j) {
             // heightNext[i][j] = ...
             double hx2=pow(x[1]-x[0],2.0);
             double hy2=pow(y[1]-y[0],2.0);
             heightNext[i][j] = dt*dt*(((u2(x[i+1],y[j])-u2(x[i-1],y[j]))*(height[i+1][j]-height[i-1][j]))/(4.0*hx2)
             + ((u2(x[i],y[j+1])-u2(x[i],y[j-1]))*(height[i][j+1]-height[i][j-1]))/(4.0*hy2)+ 
             u2(x[i],y[j])*((height[i+1][j]-2*height[i][j]+height[i-1][j])/(hx2)+(height[i][j+1]-2.0*height[i][j]+height[i][j-1])/(hy2))) 
             +2.0*height[i][j]-heightPrev[i][j];
        }
    }
}

int
main(int argc, char* argv[])
{
    // USAGE: name-of-binary [configuration-file] [<settings-to-overwrite> ...]

    // Read the default input
    string inputPath = "configuration.in";
    // Optionally override configuration file.
    if (argc > 1) inputPath = argv[1];

    ConfigFile configFile(inputPath);
    // Override settings
    for (int i = 2; i < argc; i++)
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Read inputs
    const double Lx = configFile.get<double>("Lx");
    const double Ly = configFile.get<double>("Ly");

    BoundaryConditions boundaryConditions;
    boundaryConditions.left = boundaryConditionFromString(configFile.get<string>("BoundaryL"));
    boundaryConditions.up = boundaryConditionFromString(configFile.get<string>("BoundaryU"));
    boundaryConditions.down = boundaryConditionFromString(configFile.get<string>("BoundaryD"));
    boundaryConditions.right = boundaryConditionFromString(configFile.get<string>("BoundaryR"));

    const double A = configFile.get<double>("A");
    const double omega = configFile.get<double>("omega");

    InitialCondition initialCondition =
      initialConditionFromString(configFile.get<string>("InitialCondition"));
    const double f0 = configFile.get<double>("f0");
    const double xa = configFile.get<double>("xa");
    const double xb = configFile.get<double>("xb");
    const double ya = configFile.get<double>("ya");
    const double yb = configFile.get<double>("yb");
    const double eigenmodeAmplitude = configFile.get<double>("eigenmodeAmplitude");
    const int eigenmodeM = configFile.get<int>("eigenmodeM");
    const int eigenmodeN = configFile.get<int>("eigenmodeN");

    const std::string propagationFunctionName = configFile.get<string>("PropagationFunction");
    const double u0 = configFile.get<double>("u0");
    const double g = configFile.get<double>("g");
    const double h0 = configFile.get<double>("h0");
    const double h1 = configFile.get<double>("h1");
    const double a = configFile.get<double>("a");
    const double b = configFile.get<double>("b");

    const int Nx = configFile.get<int>("Nx");
    const int Ny = configFile.get<int>("Ny");
    const double betaCFL = configFile.get<double>("betaCFL");
    const double endTime = configFile.get<double>("endTime");
    const string outputFilePath = configFile.get<string>("output");

    // TODO: Implement finite difference grid
    vector<double> x(Nx);
    vector<double> y(Ny);
    for (unsigned int i = 0; i < x.size(); ++i) x[i] = i*Lx/(x.size()-1);
    for (unsigned int j = 0; j < y.size(); ++j) y[j] = j*Ly/(y.size()-1);

    // The 'u2FunctionGenerator' returns a function that can depend on the different initial
    // conditions. After this point you can use 'u2' as a normal function, e.g. `u2(x[i], y[j])`
    // but it also has the benefit of being able to be passed as an argument to another
    // function, see e.g. 'advanceHeightField'
    WaveSpeedFunction u2 = u2FunctionGenerator(propagationFunctionName, u0, g, h0, h1, a, b, Ly);

    const double dt = calculateDt(betaCFL, x, y, u2);

    vector<vector<double>> height(Nx, vector<double>(Ny));
    vector<vector<double>> heightPrev(Nx, vector<double>(Ny));
    //cout<<"DEBUG:: avant Initial condition"<<endl;
    setInitialConditions(height,
                         heightPrev,
                         initialCondition,
                         x,
                         y,
                         dt,
                         u2,
                         f0,
                         xa,
                         xb,
                         ya,
                         yb,
                         Lx,
                         Ly,
                         eigenmodeAmplitude,
                         eigenmodeM,
                         eigenmodeN);
    vector<vector<double>> heightNext(Nx, vector<double>(Ny));

	//cout<<"DEBUG:: Finit Initial condition"<<endl;
    // Temporal loop
    const int frameCount = endTime / dt;
    cout << "2D WAVE SIMULATOR | " << inputPath << " | dt=" << setprecision(2) << dt
         << " | frames=" << frameCount << endl;
    int printCounter = 0;
    //cout<<"DEBUG:: avant d'entre dans la boucle"<<endl;
    
    //writeData(getFilepathForFrame(0, outputFilePath), 0.0, x, y, height);
    for (int frame = 1; frame <= frameCount; ++frame) {
        double time = frame * dt;
		
		//cout<<"DEBUG:: avant advanceHeightField"<<endl;
        advanceHeightField(time, dt, heightNext, height, heightPrev, x, y, u2);
        ///cout<<"DEBUG:: après advanceHeightField"<<endl;
        setBoundaryConditions(boundaryConditions, time, dt, heightNext, height, x, y, A, omega, u2);
		///cout<<"DEBUG:: apres setBoundaryConditions"<<endl;
		
        // Rotate the buffers
        std::swap(heightPrev, height);
        std::swap(height, heightNext);

        const string framepath = getFilepathForFrame(frame, outputFilePath);
        if (verbose) {
            // Print every 10%
            const double completedPercentage = std::round(time / endTime * 100.0);
            if (int(completedPercentage) / 10 > printCounter) {
                printCounter = int(completedPercentage) / 10;
                std::cout << "  Frame " << to_string(frame) << " (" << setprecision(3)
                          << completedPercentage << "%) | dt=" << dt << " | output=" << framepath
                          << std::endl;
            }
        }
		//cout<<"DEBUG:: avant writeData"<<endl;
        writeData(framepath, time, x, y, height);
        //cout<<"DEBUG:: apres writeData"<<endl;
    }
	//cout<<"DEBUG:: main return 0"<<endl;
    return 0;
}

