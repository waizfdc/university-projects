//проверить, была ли стрелка в пути
//поиск красной стрелки на хардах
//врубить коррекцию
//рисовать прямоугольник у  клада

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <vector>
#include <tuple>

using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::ofstream;
using std::vector;
using std::tuple;
using std::tie;
using std::make_tuple;

#include "io.h"
#include "matrix.h"

enum object_type {none = 0, arrow = 1, treasure = 2};

//порог бинаризации. возможно, нужно вычислять вручную
const uint binarisation_threshold = 150;

//нахер не нужно
//typedef tuple<uint, uint, uint, uint> Rect;

//картиночка для вывода. пока для отладки использую
Image to_print;

//поток вывода в файл. нужно задавать вручную
//ofstream fout;

class BoxFilterOp
{
public:
    BoxFilterOp() : radius(2) {}
    BoxFilterOp(int rad) : radius(rad) {}
    tuple<uint, uint, uint> operator () (const Image &m) const
    {
        uint size = 2 * radius + 1;
        uint r, g, b, sum_r = 0, sum_g = 0, sum_b = 0;
        for (uint i = 0; i < size; ++i) {
            for (uint j = 0; j < size; ++j) {
                // Tie is useful for taking elements from tuple
                tie(r, g, b) = m(i, j);
                sum_r += r;
                sum_g += g;
                sum_b += b;
            }
        }
        auto norm = size * size;
        sum_r /= norm;
        sum_g /= norm;
        sum_b /= norm;
        return make_tuple(sum_r, sum_g, sum_b);
    }
    // Radius of neighbourhoud, which is passed to that operator
    int radius;
};

//расстояние между точками
inline double length(uint begin_x, uint begin_y, uint end_x, uint end_y) {
    return pow(abs(end_x - begin_x), 2) + pow(abs(end_y - begin_y), 2);
}

//алгоритм брезенхэма для рисования линии
void draw_line(Image& img, uint begin_x, uint begin_y, uint end_x, uint end_y) {
    
    const bool step = (abs(end_y - begin_y) > abs(end_x - begin_x));
    if (step) {
        std::swap(begin_x, begin_y);
        std::swap(end_x, end_y);
    }

    if (begin_x > end_x) {
        std::swap(begin_x, end_x);
        std::swap(begin_y, end_y);
    }

    const float dx = end_x - begin_x;
    const float dy = abs(end_y - begin_y);

    float error = dx / 2.0f;
    const int ystep = (begin_y < end_y) ? 1 : -1;
    uint y = begin_y;

    const uint maxX = end_x;

    for (uint x = begin_x; x < maxX; x++) {
        if(step)
            img(y, x) = make_tuple(0, 255, 0);
        else
            img(x, y) = make_tuple(0, 255, 0);

        error -= dy;
        if (error < 0) {
            y += ystep;
            error += dx;
        }
    } 

}

//класс для определения объектов на картинке
class object {
public:
    uint min_x, min_y, max_x, max_y, masscentr_x, masscentr_y, distant_x, distant_y;
    int color;
    object_type type;
    object(int col = 1):
        min_x{3000},
        min_y{3000},
        max_x{0},
        max_y{0},
        masscentr_x{0},
        masscentr_y{0},
        distant_x{masscentr_x},
        distant_y{masscentr_y},
        color{col},
        type{none}
        {}
    //записывает минимум и максимум, если они поменялись. надо поменять название
    void check_border(const uint x, const uint y) {
        min_x = min_x < x ? min_x : x;
        min_y = min_y < y ? min_y : y;
        max_x = max_x > x ? max_x : x;
        max_y = max_y > y ? max_y : y;
        return;
    }

    //посути, вычисляет границы прямоугольника, описывающего фигуру
    void initialise_borders(Matrix<int>& obj_matr) { //передаём по ссылке, тк нужно менять матрицу
        min_x = obj_matr.n_rows;
        min_y = obj_matr.n_cols;
        max_x = 0;
        max_y = 0;
        for (uint i = 0; i < obj_matr.n_rows; i++)
            for (uint j = 0; j < obj_matr.n_cols; j++) 
                if (obj_matr(i, j) == color) {
                    check_border(i, j);
                    obj_matr(i, j) = 0;
                }
        return;
    }

    tuple<uint, uint> find_distant(Matrix<int>& obj_matr, uint b_x, uint b_y) {
        uint x{b_x}, y{b_y};
        double max_length{0.0};
        for (uint i = min_x; i <= max_x; i++)
            for (uint j = min_y; j <= max_y; j++)
                if (obj_matr(i, j) == color)
                    if (max_length < length(b_x, b_y, i, j)) {
                        max_length = length(b_x, b_y, i, j);
                        x = i;
                        y = j;
                    }
        return make_tuple(x, y);
    }

    //делает подготовительную работу и определяет тип объекта
    void name_type(Matrix<int>& obj_matr) {
        //проверять только точки из этого объекта. вроде сейчас путается
        masscentr_x = 0; masscentr_y = 0;
        int count{0};
        for (uint i = min_x; i <= max_x; i++)
            for (uint j = min_y; j <= max_y; j++)
                if (obj_matr(i, j) == color) {
                    masscentr_x += i;
                    masscentr_y += j;
                    count++;
                }
        masscentr_x /= count;
        masscentr_y /= count;
        //выводим центр масс
        // for (uint i = masscentr_x - 1; i < masscentr_x + 2; i++) 
        //     for  (uint j = masscentr_y -1; j < masscentr_y + 2; j++)
        //         to_print(i, j) = make_tuple(0, 0, 255);

        //поиск самой удалённой точки. вынести в отдельную функцию
        tie(distant_x, distant_y) = find_distant(obj_matr, masscentr_x, masscentr_y);
        //cout << "max_length = " << max_length << endl;
        //поиск самой удалённой точки от самой удалённой точки
        uint second_distant_x, second_distant_y;
        tie(second_distant_x, second_distant_y) = find_distant(obj_matr, distant_x, distant_y);
                
        //cout << max_length << endl;
        double first_param = 0, second_param = 0;
        for (uint i = min_x; i <= max_x; i++)
            for (uint j = min_y; j <= max_y; j++)
                if (obj_matr(i, j)) {
                    first_param += pow (length(distant_x, distant_y, i, j), 3);
                    second_param += pow (length(second_distant_x, second_distant_y, i, j), 3);
                }

        // for (uint i = second_distant_x - 5; i < second_distant_x + 6; i++) 
        //     for  (uint j = second_distant_y - 5; j < second_distant_y + 6; j++)
        //         to_print(i, j) = make_tuple(255, 255, 0);
        
        // cout << "distant 1: (" << distant_x << ", " << distant_y << ") distant 2: (" <<
        //     second_distant_x << ", " << second_distant_y << ") 1 param = " <<
        //     first_param << "; 2 param = " << second_param << endl;
        if (second_param > first_param) {
            distant_x = second_distant_x;
            distant_y = second_distant_y;
        }
        //выводим удалённую точку 
        // for (uint i = distant_x - 1; i < distant_x + 2; i++) 
        //     for  (uint j = distant_y -1; j < distant_y + 2; j++)
        //         to_print(i, j) = make_tuple(255, 0, 255);

        //draw_line(to_print, masscentr_x, masscentr_y, distant_x, distant_y);

        int x = distant_x - masscentr_x;
        x /= 2;
        int y = distant_y - masscentr_y;
        y /= 2;

        float alfa = 0.5;
        //std::cout << masscentr_x << " " << masscentr_x - x << " " << masscentr_y << " " << masscentr_y + y << " " << distant_y << endl;
        // to_print(masscentr_x - y + alfa * x, masscentr_y + x + alfa * y) = make_tuple(0, 0, 255);
        // to_print(masscentr_x + y + alfa * x, masscentr_y - x + alfa * y) = make_tuple(0, 0, 255);
        // to_print(masscentr_x - y - alfa * x, masscentr_y + x - alfa * y) = make_tuple(0, 0, 255);
        // to_print(masscentr_x + y - alfa * x, masscentr_y - x - alfa * y) = make_tuple(0, 0, 255);
        
        //проверка по четырём точкам. возможно нужно заменить
        if (obj_matr(masscentr_x - y + alfa * x, masscentr_y + x + alfa * y) && obj_matr(masscentr_x + y + alfa * x, masscentr_y - x + alfa * y) &&
            !obj_matr(masscentr_x - y - alfa * x, masscentr_y + x - alfa * y) && !obj_matr(masscentr_x + y - alfa * x, masscentr_y - x - alfa * y))
            type = arrow;
        else
            type = treasure;
        // if (type == arrow) 
        //     for (uint i = masscentr_x - 5; i < masscentr_x + 6; i++) 
        //         for  (uint j = masscentr_y - 5; j < masscentr_y + 6; j++)
        //             to_print(i, j) = make_tuple(255, 255, 0);
    }

    //возвращает true, если объект меньше 40 пикселей. зануляет его. возможно, нужно вручную задавать порог пикселей
    bool del_lit(Matrix<int>& obj_matr) {
        Matrix<int> tmp = obj_matr.deep_copy();
        int count{0};
        for (uint x = min_x; x < max_x; x++)
            for (uint y = min_y; y < max_y; y++)
                if (obj_matr(x, y) == color) {
                    tmp(x, y) = 0;
                    count ++;
                }
        if (count < 40) {
            obj_matr = tmp.deep_copy();
            return true;
        }
        return false;
    }
};

//выдаёт яркость пикселя по супер формуле. написать альтернативный вариант с максимальной яркостью по каналам. и сменить название
uint brightness_f(std::tuple<uint, uint, uint> c);
uint brightness_c(std::tuple<uint, uint, uint> c) {
    uint r, g,  b;
    tie(r, g, b) = c;
    uint max = r;
    if (g > max) max = g;
    if (b > max) max = b;
    return max;
}

//печатает значения каналов. нахуй не нужна
// void show(const Image& img) {
//     for (uint i = 0; i < img.n_rows; ++i) {
//         for (uint j = 0; j < img.n_cols; ++j) {
//             uint a, b, c;
//             tie(a, b, c) = img(i, j);
//             fout << '(' << a << '.' << b << '.' << c << ")  ";
//         }
//         fout << endl;
//     }
// }

//возвращает бинаризованное изображение
Image binarise(const Image& in) {
    Image binar = in.deep_copy();
    for (uint i = 0; i < binar.n_rows; ++i)
        for (uint j = 0; j < binar.n_cols; ++j) {
            uint r, g, b;
            tie(r, g, b) = binar(i, j);
            //поменять на функцию яркости
            r = g = b = (brightness_c(binar(i, j)) > binarisation_threshold) ? 255 : 0;
            binar(i, j) = make_tuple(r, g, b);
        }
    return binar;
}

//типа гамма-коррекция
Image gamma_cor(const Image& in) {
    double ga = 1.7;
    double c = 255 / pow(255, ga);
    Image cor = in.deep_copy();
    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j < in.n_cols; j++) {
            uint r, g, b;
            tie(r, g, b) = in(i, j);
            uint r1 = c * pow(r, ga);
            uint g1 = c * pow(g, ga);
            uint b1 = c * pow(b, ga);
            cor(i, j) = make_tuple(r1, g1, b1);  
        }
    return cor;
}

Image log_cor(const Image& in) {
    double c = 255 / log(255);
    Image cor = in.deep_copy();
    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j < in.n_cols; j++) {
            uint r, g, b;
            tie(r, g, b) = in(i, j);
            uint r1 = c * log(r);
            uint g1 = c * log(g);
            uint b1 = c * log(b);
            cor(i, j) = make_tuple(r1, g1, b1);  
        }
    return cor;
}

Image linear_cor( const Image& in) {
    uint min = 255, max = 0;

    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j< in.n_cols; j++) {
            if (brightness_f(in(i, j)) < min) min = brightness_f(in(i, j));
            if (brightness_f(in(i, j)) > max) max = brightness_f(in(i, j));
        }

    if (min <  10) min = 13;
    if (max > 245) max = 240;
    double c = 255 / (max - min);

    Image res = in.deep_copy();

    for (uint i = 0; i < in.n_rows; i++)
        for (uint j = 0; j< in.n_cols; j++) {
            int r, g, b;
            tie(r, g, b) = in(i, j);
            r -= min * 0.2989; if (r < 0) r = 0;
            g -= min * 0.5870; if (g < 0) g = 0;
            b -= min * 0.1140; if (b < 0) b = 0;

            r *= c; if (r > 255) r = 255;
            g *= c; if (g > 255) g = 255;
            b *= c; if (b > 255) b = 255;

            uint r1 = r; uint g1 = g; uint b1 = b;

            res(i, j) = make_tuple(r1, g1, b1);
        }
        return res;
}

//tuple<vector<Rect>, Image>
// find_treasure(const Image& in)
// {
//     // Base: return Rect of treasure only
//     Image processed = in.unary_map(BoxFilterOp(2));
//     processed = binarise(processed);

//     // Bonus: return Rects of arrows and then treasure Rect
//     auto path = vector<Rect>();
//     return make_tuple(path, in.deep_copy());
// }

//поменять название
uint brightness_f(std::tuple<uint, uint, uint> c) {
    int r, g, b;
    tie(r, g, b) = c;
    double gray = 0.2989 * r + 0.5870 * g + 0.1140 * b;
    return gray;
}

//какая-то херня. вроде не нужна
// int grayworld_pr(const Image& in) {

//     double gray = 0;

//     for (uint i = 0; i < in.n_rows; ++i)
//         for (uint j = 0; j < in.n_cols; ++j) 
//             gray += brightness_f(in(i, j));

//         gray /= in.n_cols * in.n_rows;
//         return gray;
// }

//моя суперская коррекция. вычетает размытие из оригинала. работает долго. добавить параметр - радиус размытия
Image my_cor(const Image& in, uint rad) {
    Image res = in.unary_map(BoxFilterOp(rad));
    for (uint i = 0; i < res.n_rows; i++)
        for (uint j = 0; j< res.n_cols; j++) {
            int r1, g1, b1;
            uint r2, g2, b2;
            tie(r1, g1, b1) = in(i, j);
            tie(r2, g2, b2) = res(i, j);
            r1 -= r2; if (r1 < 0) r1 = 0;
            g1 -= g2; if (g1 < 0) g1 = 0;
            b1 -= b2; if (b1 < 0) b1 = 0;
            r2 = r1;  g2 = g1; b2 = b1;
            res(i, j) = make_tuple(r2, g2, b2);
        }
        return res;
    
}

//размечает бинаризованное изображение, выделяет компоненты связности. работает. супер алгоритм с хабра
Matrix<int> find_objects(const Image& in) {
    auto field = Matrix<int>(in.n_rows, in.n_cols);
    // std::cout << "in funk" << endl;
    for (uint i = 0; i < in.n_rows; i++) 
        for (uint j = 0; j < in.n_cols; j++) {
            uint r, g, b;
            tie(r, g, b) = in(i, j);
            field(i, j) = r ? 1 : 0;
        }

    int teg = 2;

    if (field(0, 0)) 
        field(0, 0) = teg++;
        // std::cout << "2" << endl;
    for (uint i = 1; i < field.n_cols; i++)

        if (field(0, i)) 
            field(0, i) = (field(0, i - 1)) ? field(0, i - 1) : teg++;

    for (uint i = 1; i < field.n_rows; i++) {

        if (field(i, 0))
            field(i, 0) = (field(i - 1, 0)) ? field(i - 1, 0) : teg++;

        for (uint j = 1; j < field.n_cols; j++) {
            int a = field(i, j), b = field(i, j - 1), c = field(i - 1, j);
            if (a) {
                if (!b && !c)
                    field(i, j) = teg++;
                if (!b && c)
                    field(i, j) = c;
                if (b && !c) 
                    field(i, j) = b;
                if (b && c) {
                    field (i, j) = c;
                    if (b != c) 
                        for (uint k = 0; k <= i; k++)
                            for (uint l = 0; l < field.n_cols; l++)
                                if (field(k, l) == b) 
                                    field(k, l) = c;
                }
            }
        }
    }

    return field;
}

//ща вроде нахуй не нужны
// std::vector<object> arrows;
// std::vector<object> tresh;
// std::vector<object> all;

//закинывает объекты в вектор. возвращать этот вектор, вместо заполнения конкретного
std::vector<object> get_objects(const Image& binarised_im) {
    std::vector<object> all;
    Matrix<int> obj_matr = find_objects(binarised_im);
    Matrix<int> obj_field = obj_matr.deep_copy();
    for (uint i = 0; i < obj_matr.n_rows; i++) 
        for (uint j = 0; j < obj_matr.n_cols; j++) 
            if(obj_matr(i, j)) {
                object new_obj = object(obj_matr(i, j));
                new_obj.initialise_borders(obj_matr);
                if (new_obj.del_lit(obj_field))
                    continue;
                new_obj.name_type(obj_field);
                all.push_back(new_obj);
                // if(new_obj.type == arrow)
                //     arrows.push_back(new_obj);
                // else
                //     tresh.push_back(new_obj);
            }
    return all;
}

//находит и рисует путь
void make_path(const Image& img, std::vector<object> all) {
    object current, next_obj;
    //допилить поиск красной стрелки, чтоб работал на харде
    {
        uint r, g, b;
        for (const auto &obj : all) {
            tie(r, g, b) = img(obj.masscentr_x, obj.masscentr_y);
            if (obj.type == arrow && r > 225 && g < 30 && b < 30) {
                current = obj;
                break;
            }
        }
    }
    //напечатать текущий объект в файл

    //составление пути
    while (current.type == arrow) {
        //идём с шагом в четверь длины стрелки и смотрим, в объекте ли мы
        int vect_x = current.distant_x - current.masscentr_x, vect_y = current.distant_y - current.masscentr_y;
        vect_x /= 2; vect_y /= 2;
        //добавить проверку наа отрицательность cur_x cur_y
        //если уж они станут отрицательными, то станут очень большими => не пройдут проверку на величину. так что не надо
        uint cur_x = current.distant_x + vect_x, cur_y = current.distant_y + vect_y;
        if (cur_x > img.n_rows || cur_y > img.n_cols) {
            std::cout << "error " << vect_x << " " << vect_y<< endl;
            return;
        }
        bool found = false;
        //проверка на отрицательность. аналогично не нужна
        while (!(cur_x > img.n_rows || cur_y > img.n_cols)) {
           
            for (const auto &obj : all) {
                if (cur_x > obj.min_x && cur_y > obj.min_y && cur_x < obj.max_x && cur_y < obj.max_y) {
                    next_obj = obj;
                    found = true;
                    break;
                }
            }
            if (found) {
                draw_line(to_print, current.masscentr_x, current.masscentr_y, next_obj.masscentr_x, next_obj.masscentr_y);
                //записать координаты следующего в файл
                current = next_obj;
                break;
            }
            //+проверка на отрицательность. возможно сначала приравнивать к int переменной и сравнивать ёё с нулём
            cur_x += vect_x;
            cur_y += vect_y;
        }
        if (!found) {
            std::cout << "error" << endl;
            return;
        }

    }
    //рисуем прямоугольник вокруг current
}

int main(int argc, char **argv)
{
    // fout << "hyi" << endl;
    // fout << "pizda" << endl;
    ofstream fout(argv[3]);
    if (!argc) return 0;

    Image src_image = load_image(argv[1]);
    to_print = src_image.deep_copy();
    Image binar = linear_cor(src_image);
    binar = gamma_cor(binar);
    binar = binar.unary_map(BoxFilterOp(2));
    binar = binarise(binar);
    auto all = get_objects(binar);
    make_path(src_image, all);
    // for (uint i = 0; i < objs.size(); i++){
    //     for (uint k = objs[i].min_x; k < objs[i].max_x; k++) {
    //         src_image(k, objs[i].min_y) = make_tuple(0, 255, 0);
    //         src_image(k, objs[i].max_y) = make_tuple(0, 255, 0);
    //     }
    //     for (uint k = objs[i].min_y; k < objs[i].max_y; k++) {
    //         src_image(objs[i].min_x, k) = make_tuple(0, 255, 0);
    //         src_image(objs[i].max_x, k) = make_tuple(0, 255, 0);
    //     }
    // }
    save_image(to_print, "output/objects.bmp");
    if (1) return 0;

    //std::cout << "ok" << endl;
    //Matrix<int> finder = find_objects(src_image);
    //std::cout << "after funk" << endl;
    //std::cout << finder << endl;
    

    // Image filt_my = my_cor(src_image, 6);
    // save_image(src_image, "output/source.bmp");
    // save_image(filt_my, "output/mycor.bmp");
    // Image filt_lin = linear_cor(src_image);
    // save_image(filt_lin, "output/linear.bmp");
    // Image filt_gamma = gamma_cor(src_image);
    // save_image(filt_gamma, "output/gamma.bmp");
    // Image filt_log = log_cor(src_image);
    // save_image(filt_log, "output/log.bmp");
    // //cout << "grayworld: source - " << grayworld_pr(src_image) << " linear - " <<
    // //grayworld_pr(filt_lin) << " gamma - " << grayworld_pr(filt_gamma) <<
    // " log - " << grayworld_pr(filt_log) << endl;
    // Image filt = gamma_cor(filt_lin);
    // Image processed = filt.unary_map(BoxFilterOp(2));
    // save_image(processed, "output/middle.bmp");
    // processed = binarise(processed);
    // save_image(processed, argv[2]);
    // filt = filt_gamma.unary_map(BoxFilterOp(14));
    // filt = binarise(filt);
    // save_image(filt, "output/gaus.bmp");

    return 0;    
    /*
    if (argc != 4)
    {
        cout << "Usage: " << endl << argv[0]
             << " <in_image.bmp> <out_image.bmp> <out_path.txt>" << endl;
        return 0;
    }

    try {
        Image src_image = load_image(argv[1]);
        ofstream fout(argv[3]);

        vector<Rect> path;
        Image dst_image;
        tie(path, dst_image) = find_treasure(src_image);
        save_image(dst_image, argv[2]);

        uint x, y, width, height;
        for (const auto &obj : path)
        {
            tie(x, y, width, height) = obj;
            fout << x << " " << y << " " << width << " " << height << endl;
        }

    } catch (const string &s) {
        cerr << "Error: " << s << endl;
        return 1;
    }
    */
}
