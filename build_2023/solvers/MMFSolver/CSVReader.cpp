#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        return 1;
    }

    char *filename = argv[1];
    ofstream csv_file;
    csv_file.open(filename, ios::in);

    vector<vector<double>> csv;

    if (csv_file.is_open())
    {
        char *row;

        while (getline(csv_file, row))
        {
            vector<double> row_vec;
            csv.push_back(row_vec);
            char *col;
            double col_d;
            col = strtok(row, ",");
            while (col != NULL)
            {
                col_d = atof(col);
                row_vec.push_back(col_d);
                col = strtok(NULL, ",");
            }
        }
    }
    else
    {
        return 1;
    }

    return 0;
}