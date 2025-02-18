#include <bits/stdc++.h>
#include "bitmap_image.hpp"

using namespace std;

typedef struct {
    double x, y, z;
    double n=1;
} Point;

typedef struct {
    double r, g, b;
} Color;

typedef Point Vector;

const Vector X_UNIT_VEC = {.x=1, .y=0, .z=0};
const Vector Y_UNIT_VEC = {.x=0, .y=1, .z=0};
const Vector Z_UNIT_VEC = {.x=0, .y=0, .z=1};

typedef struct {
    set<Point> points;
} Triangle;

double deg2rad(double deg) {
    return deg * M_PI / 180;
}

double rad2deg(double rad) {
    return rad * 180 / M_PI;
}

stack<double**> S;

int screen_width, screen_height;
double dx, dy;


double ** mult_mat(double** A, double** B) {
    double ** C = new double*[4];
    for (int i = 0; i < 4; i++) {
        C[i] = new double[4];
        for (int j = 0; j < 4; j++) {
            C[i][j] = 0;
            for (int k = 0; k < 4; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}


void print_mat(double** arg) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << arg[i][j] << " ";
        }
        cout << endl;
    }
}


void dealloc_mat(double** arg) {
    for (int i = 0; i < 4; i++) {
        delete[] arg[i];
    }
    delete[] arg;
}


double** alloc_mat(double diag = 0) {
    double** M = new double*[4];
    for (int i = 0; i < 4; i++) {
        M[i] = new double[4];
        for (int j = 0; j < 4; j++) {
            if (diag != 0 && i == j)
                M[i][j] = diag;
            else M[i][j] = 0;

        }
    }
    M[3][3] = 1;
    return M;
}


Point mult_mat_point(double** M, Point p) {
    Point res;
    res.x = M[0][0] * p.x + M[0][1] * p.y + M[0][2] * p.z + M[0][3];
    res.y = M[1][0] * p.x + M[1][1] * p.y + M[1][2] * p.z + M[1][3];
    res.z = M[2][0] * p.x + M[2][1] * p.y + M[2][2] * p.z + M[2][3];
    res.n = M[3][0] * p.x + M[3][1] * p.y + M[3][2] * p.z + M[3][3];

    res.x /= res.n;
    res.y /= res.n;
    res.z /= res.n;
    res.n /= res.n;

    return res;
}


Vector scale_to_r(Vector v, double r) {
    double len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    v.x *= r / len;
    v.y *= r / len;
    v.z *= r / len;
    return v;
}


Vector scalar_mult_vec(double r, Vector v) {
    v.x *= r;
    v.y *= r;
    v.z *= r;
    return v;
}


double dot_vec(Vector a, Vector b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}


Vector cross_vec(Vector a, Vector b) {
    Vector res;
    res.x = a.y * b.z - a.z * b.y;
    res.y = a.z * b.x - a.x * b.z;
    res.z = a.x * b.y - a.y * b.x;
    return res;
}


Vector add_vec(Vector a, Vector b) {
    Vector res;
    res.x = a.x + b.x;
    res.y = a.y + b.y;
    res.z = a.z + b.z;
    return res;
}


Vector rodrig(Vector arg, Vector axis, double theta) {
    Vector res;
    res = scalar_mult_vec(cos(deg2rad(theta)), arg);
    res = add_vec(res, scalar_mult_vec((1-cos(deg2rad(theta)))*dot_vec(axis, arg), axis));
    res = add_vec(res, scalar_mult_vec(sin(deg2rad(theta)) , cross_vec(axis, arg)));
    return res;
}


static unsigned long int g_seed = 1;

int custom_rand_int() {
    g_seed = (214013 * g_seed + 2531011);
    return (g_seed >> 16) & 0x7FFF;
}


int find_min_x_index(vector<Point> p) {
    int min_x = 0;
    // cout<<"from find min x"<<p.size()<<endl;
    int n = p.size();
    for (int i = 0; i < n; i++) {
        if (p[i].x < p[min_x].x) {
            min_x = i;
        }
    }
    return min_x;
}

int find_max_x_index(vector<Point> p) {
    int max_x = 0;
    // cout<<"from find max x"<<p.size()<<endl;
    int n = p.size();
    for (int i = 0; i < n; i++) {
        if (p[i].x > p[max_x].x) {
            max_x = i;
        }
    }
    return max_x;
}

int find_min_y_index(vector<Point> p) {
    int min_y = 0;
    int n = p.size();
    for (int i = 0; i < n; i++) {
        if (p[i].y < p[min_y].y) {
            min_y = i;
        }
    }
    return min_y;
}

int find_max_y_index(vector<Point> p) {
    int max_y = 0;
    int n = p.size();
    for (int i = 0; i < n; i++) {
        if (p[i].y > p[max_y].y) {
            max_y = i;
        }
    }
    return max_y;
}

int point_to_row(Point p) {
    double top_y = 1 - dy / 2;
    return round((top_y - p.y) / dy);
}

int point_to_col(Point p) {
    double left_x = -1 + dx / 2;
    return round((p.x - left_x) / dx);
}



Point eye, look;
Vector up;
double fov_y, aspect_ratio, near, far;

int total_triangles = 0;

int main() {
    double** M = alloc_mat(1);
    

    ifstream scene("scene.txt");
    ofstream stage1_out("stage1.txt");
    string line;


    scene >> eye.x >> eye.y >> eye.z;
    scene >> look.x >> look.y >> look.z;
    scene >> up.x >> up.y >> up.z;


    scene >> fov_y >> aspect_ratio >> near >> far;

    while (true) {
        scene>>line;
        if (line == "triangle") {
            total_triangles++;
            for (int i = 0; i < 3; i++) {
                Point p;
                scene >> p.x >> p.y >> p.z;
                p = mult_mat_point(M, p);
                stage1_out << fixed << setprecision(7);
                stage1_out << p.x << " " << p.y << " " << p.z << endl;
            }
            stage1_out << endl;

        } else if (line == "translate") {
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            double** T = alloc_mat(1);
            T[0][3] = tx;
            T[1][3] = ty;
            T[2][3] = tz;
            double** old_M = M;
            M = mult_mat(M, T);
            dealloc_mat(T);
            dealloc_mat(old_M);

        } else if (line == "scale") {
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            double** S = alloc_mat(1);
            
            S[0][0] = sx;
            S[1][1] = sy;
            S[2][2] = sz;
            double** old_M = M;
            M = mult_mat(M, S);
            dealloc_mat(S);
            dealloc_mat(old_M);

        } else if (line == "rotate") {
            double angle, rx, ry, rz;
            scene >> angle >> rx >> ry >> rz;
            double** R = alloc_mat(0);

            Vector axis = {.x=rx, .y=ry, .z=rz};
            axis = scale_to_r(axis, 1);
            Vector c1 = rodrig(X_UNIT_VEC, axis, angle);
            Vector c2 = rodrig(Y_UNIT_VEC, axis, angle);
            Vector c3 = rodrig(Z_UNIT_VEC, axis, angle);


            R[0][0] = c1.x; R[0][1] = c2.x; R[0][2] = c3.x;
            R[1][0] = c1.y; R[1][1] = c2.y; R[1][2] = c3.y;
            R[2][0] = c1.z; R[2][1] = c2.z; R[2][2] = c3.z;
            
            double** old_M = M;
            M = mult_mat(M, R);
            dealloc_mat(old_M);
            dealloc_mat(R);

        } else if (line == "push") {
            double ** M_to_push = alloc_mat(0);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4;j++) {
                    M_to_push[i][j] = M[i][j];
                }
            }
            S.push(M_to_push);

        } else if (line == "pop") {
            dealloc_mat(M);
            M = S.top();
            S.pop();

        } else if (line == "end") {
            break;

        } else {
            cout << "Invalid command: " << line << endl;

        }
    }

    dealloc_mat(M);
    stage1_out.close();
    scene.close();
    // stage 2
    
    ifstream stage1_in("stage1.txt");
    ofstream stage2_out("stage2.txt");

    Vector l_vec = {.x=look.x-eye.x, .y=look.y-eye.y, .z=look.z-eye.z};
    l_vec = scale_to_r(l_vec, 1);
    Vector r_vec = cross_vec(l_vec, up);
    r_vec = scale_to_r(r_vec, 1);
    Vector u_vec = cross_vec(r_vec, l_vec);
    u_vec = scale_to_r(u_vec, 1);

    double ** T = alloc_mat(1);
    T[0][3] = -eye.x;
    T[1][3] = -eye.y;
    T[2][3] = -eye.z;

    double ** R = alloc_mat(1);
    R[0][0] = r_vec.x;  R[0][1] = r_vec.y;  R[0][2] = r_vec.z;
    R[1][0] = u_vec.x;  R[1][1] = u_vec.y;  R[1][2] = u_vec.z;
    R[2][0] = -l_vec.x; R[2][1] = -l_vec.y; R[2][2] = -l_vec.z;


    for (int i = 0; i < total_triangles; i++) {
        // Triangle t;
        for (int j = 0; j < 3; j++) {
            Point p;
            stage1_in >> p.x >> p.y >> p.z;
            // cout<<"before T"<<endl;
            // cout << p.x << " " << p.y << " " << p.z << endl;
            p = mult_mat_point(T, p);
            // cout<<"after T"<<endl;
            // print_mat(T);
            // cout << p.x << " " << p.y << " " << p.z << endl;
            p = mult_mat_point(R, p);
            // cout<<"after R"<<endl;
            // cout << p.x << " " << p.y << " " << p.z << endl;
            stage2_out << fixed << setprecision(7);
            stage2_out << p.x << " " << p.y << " " << p.z << endl;
            // cout<<endl;
        }
        // stage1_in >> line;
        stage2_out << endl;
        // cout<<endl<<endl<<endl;

        
    }

    dealloc_mat(T);
    dealloc_mat(R);
    stage1_in.close();
    stage2_out.close();

    // stage 3

    ifstream stage2_in("stage2.txt");
    ofstream stage3_out("stage3.txt");

    double fov_x = fov_y * aspect_ratio;
    double t = near * tan(deg2rad(fov_y / 2));
    double r = near * tan(deg2rad(fov_x / 2));

    double ** P = alloc_mat();
    P[0][0] = near / r; P[0][1] = 0;        P[0][2] = 0;                        P[0][3] = 0;
    P[1][0] = 0;        P[1][1] = near/t;   P[1][2] = 0;                        P[1][3] = 0;
    P[2][0] = 0;        P[2][1] = 0;        P[2][2] = -(far+near)/(far-near);   P[2][3] = -(2*far*near)/(far-near);
    P[3][0] = 0;        P[3][1] = 0;        P[3][2] = -1;                       P[3][3] = 0;

    // print_mat(P);



    for (int i = 0; i < total_triangles; i++) {
        // Triangle t;
        for (int j = 0; j < 3; j++) {
            Point p;
            stage2_in >> p.x >> p.y >> p.z;
            // cout<<"before T"<<endl;
            // cout << p.x << " " << p.y << " " << p.z << endl;
            p = mult_mat_point(P, p);
            // cout<<"after T"<<endl;
            // print_mat(T);
            // cout << ">" << p.x << " " << p.y << " " << p.z << endl;
            // p = mult_mat_point(R, p);
            // cout<<"after R"<<endl;
            // cout << p.x << " " << p.y << " " << p.z << endl;
            stage3_out << fixed << setprecision(7);
            stage3_out << p.x << " " << p.y << " " << p.z << endl;
            // cout<<endl;
        }
        // stage1_in >> line;
        stage3_out << endl;
        // cout<<endl<<endl<<endl;

        
    }

    dealloc_mat(P);
    stage2_in.close();
    stage3_out.close();

    // stage 4
    ifstream stage3_in("stage3.txt");
    ifstream config_in("config.txt");
    ofstream z_buffer_out("z_buffer.txt");

    // read data
    config_in >> screen_width >> screen_height;

    Color bg_color = {.r=0, .g=0, .b=0};
    Color triangle_color[total_triangles];

    // initialize z-buffer
    double ** z_buffer = new double*[screen_height];
    int ** color_buffer = new int*[screen_height]; // holds the index of the triangle that is drawn at that pixel
    for (int i = 0; i < screen_height; i++) {
        z_buffer[i] = new double[screen_width];
        color_buffer[i] = new int[screen_width];
        for (int j = 0; j < screen_width; j++) {
            z_buffer[i][j] = 1;
            color_buffer[i][j] = -1;
        }
    }
    dx = (1.0 - (-1.0)) / (double)screen_width;
    dy = (1.0 - (-1.0)) / (double)screen_height;

    double top_y = 1.0 - dy / 2.0;
    double bottom_y = -1.0 + dy / 2.0;
    double left_x = -1.0 + dx / 2.0;
    double right_x = 1.0 - dx / 2.0;

    for (int i = 0; i < total_triangles; i++) {
        vector<Point> p(3);
        triangle_color[i].r = custom_rand_int() % 256;
        triangle_color[i].g = custom_rand_int() % 256;
        triangle_color[i].b = custom_rand_int() % 256;

        
        for (int j = 0; j < 3; j++) {
            stage3_in >> p[j].x >> p[j].y >> p[j].z;
            cout<<p[j].x<<" "<<p[j].y<<" "<<p[j].z<<endl;
        }
        // stage1_in >> line;
        // stage3_out << endl;
        // cout<<endl<<endl<<endl;

        int max_y_index = find_max_y_index(p);
        int min_y_index = find_min_y_index(p);

        double max_y = p[max_y_index].y;
        double min_y = p[min_y_index].y;

        int max_x_index = find_max_x_index(p);
        int min_x_index = find_min_x_index(p);

        double max_x = p[max_x_index].x;
        double min_x = p[min_x_index].x;

        double clipped_min_y = max(min_y, bottom_y);
        double clipped_max_y = min(max_y, top_y);
        double clipped_min_x = max(min_x, left_x);
        double clipped_max_x = min(max_x, right_x);

        // assumption: y axis has screen_height pixels
        // x axis has screen_width pixels

        // find the row range, top row is 0, bottom row is screen_height - 1
        int start_row = round((top_y - clipped_max_y) / dy);
        int end_row = round((top_y - clipped_min_y) / dy);

        for(int current_row = start_row; current_row <= end_row; current_row++) {
            double scanning_y = top_y - current_row * dy;
            

            vector<Point> intersections;



            // iterate through the 3 edges of the triangle
            for (int j = 0; j < 3; j++) {
                int next_j = (j + 1) % 3;
                if (p[j].y == p[next_j].y) {
                    continue;
                }

                

                if (scanning_y >= p[j].y && scanning_y <= p[next_j].y || scanning_y <= p[j].y && scanning_y >= p[next_j].y) {
                    double intersect_x = p[j].x + (p[next_j].x - p[j].x) * (scanning_y - p[j].y) / (p[next_j].y - p[j].y);
                    double intersect_z = p[j].z + (p[next_j].z - p[j].z) * (scanning_y - p[j].y) / (p[next_j].y - p[j].y);
                    
                    Point intersect_point = {.x=intersect_x, .y=scanning_y, .z=intersect_z};
                    intersections.push_back(intersect_point);
                    // cout<<"intersect: "<<intersect_x<<" "<<scanning_y<<" "<<intersect_z<<endl;
                }
            }

            if (intersections.size() < 1) continue;

            Point start_intersection = intersections[find_min_x_index(intersections)];
            Point end_intersection = intersections[find_max_x_index(intersections)];

            int start_col = point_to_col(start_intersection);
            if(start_col < 0) start_col = 0;
            else if(start_col >= screen_width) continue;
            
            int end_col = point_to_col(end_intersection);
            if(end_col < 0) continue;
            else if(end_col >= screen_width) end_col = screen_width - 1;

            for(int current_col = start_col; current_col <= end_col; current_col++) {
                double scanning_x = left_x + current_col * dx;
                double intersect_z = start_intersection.z + (end_intersection.z - start_intersection.z) * (scanning_x - start_intersection.x) / (end_intersection.x - start_intersection.x);
                if(intersect_z > 1 || intersect_z < -1) continue;
                if (intersect_z < z_buffer[current_row][current_col]) {
                    z_buffer[current_row][current_col] = intersect_z;
                    color_buffer[current_row][current_col] = i;
                    cout<<"z update "<<current_row<<" "<<current_col<<": "<<intersect_z<<endl;
                }
            }
        }

        


        cout<<endl;

        
    }

    // write to z-buffer.txt
    for (int i = 0; i < screen_height; i++) {
        for (int j = 0; j < screen_width; j++) {
            if(color_buffer[i][j] != -1) {
                z_buffer_out << fixed << setprecision(6);
                z_buffer_out << z_buffer[i][j]<<"\t";
            }
        }
        z_buffer_out << endl;
    }

    // bitmape image generation
    bitmap_image image(screen_width, screen_height);
    for (int i = 0; i < screen_height; i++) {
        for (int j = 0; j < screen_width; j++) {
            if(color_buffer[i][j] != -1) {
                Color c = triangle_color[color_buffer[i][j]];
                image.set_pixel(j, i, c.r, c.g, c.b);
            } else {
                image.set_pixel(j, i, bg_color.r, bg_color.g, bg_color.b);
            }
        }
    }
    image.save_image("out.bmp");

    for (int i = 0; i < screen_height; i++) {
        delete[] z_buffer[i];
        delete[] color_buffer[i];
        
    }

    delete[] z_buffer;
    delete[] color_buffer;

    stage3_in.close();
    config_in.close();
    z_buffer_out.flush();
    z_buffer_out.close();


    return 0;
}
