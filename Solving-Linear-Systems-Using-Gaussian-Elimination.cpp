#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <limits>

using namespace std;

using Matrix = vector<vector<double>>;
using Vector = vector<double>;

const double EPSILON = 1e-12;

struct LinearSystemResult{
    Vector solution;
    double determinant;
    double residualNorm;
};

double calculateResidualNorm(const Matrix& A, const Vector& x, const Vector& b){
    int n = A.size();
    double maxResidual = 0.0;

    for (int i = 0; i < n; i++){
        double sum = 0.0;

        for (int j = 0; j < n; j++){
            sum += A[i][j] * x[j];
        }

        maxResidual = max(maxResidual, abs(sum - b[i]));
    }

    return maxResidual;
}

LinearSystemResult solveGaussianElimination(Matrix A, Vector b){
    int n = A.size();

    if (n == 0 || A[0].size() != n || b.size() != n){
        throw invalid_argument("Invalid matrix or vector size.");
    }

    Vector scale(n);

    for (int i = 0; i < n; i++){
        scale[i] = 0.0;

        for (int j = 0; j < n; j++){
            scale[i] = max(scale[i], abs(A[i][j]));
        }

        if (scale[i] < EPSILON){
            throw runtime_error("The matrix is singular.");
        }
    }

    double determinantSign = 1.0;

    for (int k = 0; k < n - 1; k++){
        int pivotRow = k;
        double maxRatio = 0.0;

        for (int i = k; i < n; i++){
            double ratio = abs(A[i][k]) / scale[i];

            if (ratio > maxRatio){
                maxRatio = ratio;
                pivotRow = i;
            }
        }

        if (abs(A[pivotRow][k]) < EPSILON){
            throw runtime_error("The system has no unique solution.");
        }

        if (pivotRow != k) {
            swap(A[k], A[pivotRow]);
            swap(b[k], b[pivotRow]);
            swap(scale[k], scale[pivotRow]);
            determinantSign *= -1.0;
        }

        for (int i = k + 1; i < n; i++){
            double factor = A[i][k] / A[k][k];

            A[i][k] = 0.0;

            for (int j = k + 1; j < n; j++){
                A[i][j] -= factor * A[k][j];
            }

            b[i] -= factor * b[k];
        }
    }

    if (abs(A[n - 1][n - 1]) < EPSILON){
        throw runtime_error("The system has no unique solution.");
    }

    Vector x(n);

    for (int i = n - 1; i >= 0; i--){
        double sum = b[i];

        for (int j = i + 1; j < n; j++){
            sum -= A[i][j] * x[j];
        }

        x[i] = sum / A[i][i];
    }

    double determinant = determinantSign;

    for (int i = 0; i < n; i++){
        determinant *= A[i][i];
    }

    return {x, determinant, 0.0};
}

int main(){
    try {
        int n;

        cout << "Enter the number of equations: ";
        cin >> n;

        if (n <= 0){
            throw invalid_argument("The number of equations must be positive.");
        }

        Matrix A(n, Vector(n));
        Vector b(n);

        cout << "\nEnter the coefficient matrix A:" << endl;

        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                cout << "A[" << i << "][" << j << "] = ";
                cin >> A[i][j];
            }
        }

        cout << "\nEnter the constant vector B:" << endl;

        for (int i = 0; i < n; i++){
            cout << "B[" << i << "] = ";
            cin >> b[i];
        }

        Matrix originalA = A;
        Vector originalB = b;

        LinearSystemResult result = solveGaussianElimination(A, b);

        result.residualNorm = calculateResidualNorm(originalA, result.solution, originalB);

        cout << fixed << setprecision(8);

        cout << "\nSolution:" << endl;

        for (int i = 0; i < n; i++){
            cout << "x" << i + 1 << " = " << result.solution[i] << endl;
        }

        cout << "\nDeterminant of A: " << result.determinant << endl;
        cout << "Residual infinity norm: " << result.residualNorm << endl;

    } catch(const exception& error){
        cout << "\nError: " << error.what() << endl;
        return 1;
    }
    return 0;
}
