#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

// Read and write from the text file
ifstream inputFile;
ofstream outFile;

float **x_y_coordinate, **chromosomes, *fitnes, change = 1;
int dataPoints, degree , testCase , numChrom = 4, T = 3, t = 1,firstSelChro, secondSelChro;

void allocate_all_memory();
void free_all_memory();
void reedFromFile();
void initializedChromosomes();
void mean_square_Error();
void selection();
void crossover();
void mutation();
void get_fitness();
void print();
int main() {
    //inputFile.open("D:\\Semester7\\softComputing\\Assignments\\Ass_2\\20200798_20200900\\lab.txt");
    inputFile.open("D:\\Semester7\\softComputing\\Assignments\\Ass_2\\20200798_20200900\\curve_fitting_input.txt");
    outFile.open("D:\\Semester7\\softComputing\\Assignments\\Ass_2\\20200798_20200900\\curve_fitting_output.txt");

    // reed from file
    inputFile >> testCase; cout << "\nTest Case " << testCase;
    for (int i = 0; i < testCase; ++i) {
        inputFile >> dataPoints >> degree ;
        outFile << "Test Case " << i << " degree = " << degree;
        cout << "Test Case " << i << " degree = " << degree;
        allocate_all_memory();
        reedFromFile();
        initializedChromosomes();
        t = 0;
        while(t <= T) {
            for (int j = 0; j < dataPoints; ++j) {
                cout << "\n<<generations #" << t << ">>\n";
                mean_square_Error();
                selection();
                crossover();
                mutation();
            }
            if (t > T || !change){
                get_fitness(); break;
            }
            t++;
        }
        outFile << "\n\n--------------------------------------\n\n";
       // free_all_memory();
    }
    outFile.close();
    return 0;
}

void allocate_all_memory(){
    x_y_coordinate = new float* [dataPoints];
    for (int i = 0; i < dataPoints; i++) {
        x_y_coordinate[i] = new float [2];
    }

    chromosomes = new float* [numChrom];
    for (int i = 0; i < numChrom; i++) {
        chromosomes[i] = new float [degree+1];
    }

    fitnes = new float [numChrom];
}

void reedFromFile(){
    for (int i = 0; i < dataPoints; ++i)
        for (int j = 0; j < 2; ++j)
            inputFile >> x_y_coordinate[i][j];
}

void initializedChromosomes(){
    srand(time(nullptr));
    for (int i = 0; i < numChrom; ++i) {
        for (int j = 0; j < degree+1; ++j) {
            chromosomes[i][j] = ((float)rand()/RAND_MAX) * 20 - 10;
        }
    }
}

float error_at_gene(float gene, float x_coordinate, int deg){
    return (gene * pow(x_coordinate,deg));
}

float error_at_chromosome(float chromo[], float x_coordinate, float y_coordinate){
    float sum_error_chromosome = 0.0;
    for (int i = 0; i < degree + 1; ++i) {
        sum_error_chromosome += error_at_gene(chromo[i], x_coordinate, i);
    }
    return pow((sum_error_chromosome - y_coordinate),2);
}

//void mean_square_Error(){
//    float sum = 0.0, chro[4][3] = {1.95, 8.16, -2,4.26, -7.4, -2.5,  3.36, -0.3,-6.2, 0.23, 0.12, 4.62};
//
//    for (int i = 0; i < numChrom; ++i) {
//        cout << "\nerror_at_chromosome C" << i+1 ;
//        for (int j = 0; j < dataPoints; ++j) {
//            cout << "\n(" << x_y_coordinate[j][0] << " , " << x_y_coordinate[j][1] << ") = " << error_at_chromosome(chro[i], x_y_coordinate[j][0], x_y_coordinate[j][1]);
//            sum+= error_at_chromosome(chro[i], x_y_coordinate[j][0], x_y_coordinate[j][1]);
//        } fitnes[i] = (1.0 / (sum / dataPoints)); cout << "\nTotal error = " << fitnes[i];
//    }
//}

void mean_square_Error(){
    float sum = 0.0;
    for (int i = 0; i < numChrom; ++i) {
        //cout << "\nerror_at_chromosome C" << i+1 ;
        for (int j = 0; j < dataPoints; ++j) {
            //cout << "\n(" << x_y_coordinate[j][0] << " , " << x_y_coordinate[j][1] << ") = " << error_at_chromosome(chromosomes[i], x_y_coordinate[j][0], x_y_coordinate[j][1]);
            sum+= error_at_chromosome(chromosomes[i], x_y_coordinate[j][0], x_y_coordinate[j][1]);
        } fitnes[i] = (1.0 / (sum / dataPoints));
        //cout << "\nTotal error = " << fitnes[i];
    }
}

float cumulativeFitness(int index){
    float sum = fitnes[0];
    for (int i = 0; i < index; ++i) {
        sum+=fitnes[i];
    }
    return sum;
}


int selectChromosome(float r){
    if(r >= 0 && r <= cumulativeFitness(0))
        return 0;

    for (int i = 1, j = 0; i <= numChrom - 1; i++, j++) {
        if(r >= cumulativeFitness(j) && r <= cumulativeFitness(i))
            return j;
    }
}

float getRandom(){
    int lastFitness = numChrom - 1;
    float cumLastFit = cumulativeFitness(lastFitness);
    return ((float)rand()/(RAND_MAX+1.0)) * (cumLastFit - 0)+0;
}

void selection(){
    srand(time(nullptr));
    float r1 , r2;
    r1 = getRandom();
    r2 = getRandom();
    firstSelChro = selectChromosome(r1);
    secondSelChro = selectChromosome(r2);
    cout << "\nr1 = " << r1 << "\nr2 = " << r2;
    while (firstSelChro == secondSelChro){r2 = getRandom(); secondSelChro = selectChromosome(r2);}
    cout << "\n<Step4> Selection \n r1 = " << r1 << " r2 = " << r2 << "\n firstSelChro = C" << firstSelChro << "\n secondSelChro = C" << selectChromosome <<endl;
}

void crossover(){
    cout << "\n<Step5> Crossover between " << chromosomes[firstSelChro][0] << "  " << chromosomes[secondSelChro][0];


    srand(time(nullptr));

    int Xc =  rand() % ((degree) - 1 + 1) + 1;
    float  r2 ,pc;
    cout << "\n r1 = " << Xc;
    r2 =  0.0 + static_cast<float>(rand() / (RAND_MAX/((float)(1.0) - 0.0)));
    cout << "\n r2 = " << r2;
    pc =  0.5 + static_cast<float>(rand() / (RAND_MAX/((float)(0.7 - 0.4) - 0.4)));
    cout << "\n pc = " << pc;

    // check constrain of crossover
    if (r2 <= pc){
        change = 0;
        swap(chromosomes[firstSelChro][Xc], chromosomes[secondSelChro][Xc]);
        cout << "\n newChromosome :  [ "; cout << chromosomes[firstSelChro][Xc] << "  " << chromosomes[secondSelChro][Xc] << "]";
    } else {
        change = 0;
        cout << "\n no crossover";
    }
}


void non_uniform_mutation(int index ,float rm, float  pm){
    srand(time(nullptr));
    float dl, du,dy, r1, y, r;
    for (int i = 0; i < degree+1; ++i) {
        dl = chromosomes[index][i] - (-10);
        du = 10 - chromosomes[index][i];
        r1 = 0.0 + static_cast<float>(rand() / (RAND_MAX/((float)(1.0 - 0.0) - 0.0)));
        if(r1 <= 0.5)
            y = dl;
        else
            y = du;

        r = 0.0 + static_cast<float>(rand() / (RAND_MAX/((float)(1.0 - 0.0) - 0.0)));
        dy = y * (1 - (pow(r, pow((1-((float)(t/T))), 3))));
        if(y == dl)
            chromosomes[index][i] -= dy;
        else
            chromosomes[index][i] = dy - chromosomes[index][i];
    }
}

void mutation(){
    srand(time(nullptr));
    float rm, pm;
    pm = 0.001 + static_cast<float>(rand() / (RAND_MAX/((float)(0.1 - 0.001) - 0.001)));
    rm = 0.0 + static_cast<float>(rand() / (RAND_MAX/((float)(1.0 - 0.0) - 0.0)));

    if(rm <= pm) {
        change = 1;
        non_uniform_mutation(firstSelChro, rm, pm);
        cout << "\n new chromosome O1 = " << chromosomes[firstSelChro][0] << "\n";
        non_uniform_mutation(secondSelChro, rm, pm);
        cout << "\n new chromosome O2 = " << chromosomes[secondSelChro][0] << "\n";
    }
    else{change = 0;
        cout << "\n no mutation";}

}

int get_index_of_mini_chromo(){
    float mini = fitnes[0]; int i,c;
    for (i = 0, c = 0; i < numChrom; ++i) {
        if(fitnes[i] < mini){
            mini = fitnes[i];
            c = i;
        }
    }
    return c;
}

void get_fitness(){
    int fit_chromo = get_index_of_mini_chromo();
    outFile << "\ncoefficients: [ ";
    for (int i = 0; i < degree + 1; ++i) {
        outFile << chromosomes[fit_chromo][i] << " ";
    }
    outFile << "]\n cost: " << fitnes[fit_chromo];
}


void free_all_memory(){
    for (int i = 0; i < dataPoints; ++i)
        delete []x_y_coordinate[i];

    for (int i = 0; i < degree + 1; ++i)
        delete []chromosomes[i];

    delete []chromosomes;
    delete []x_y_coordinate;
    delete[]fitnes;
}



void print(){
    for (int i = 0; i < numChrom; ++i) {
        for (int j = 0; j < degree + 1; ++j) {
            cout << chromosomes[i][j] << " ";
        }cout << endl;
    }
}