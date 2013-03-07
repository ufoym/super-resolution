#pragma once
#include <fstream>

//-------------------------------------------------------------------------
// image class
//-------------------------------------------------------------------------

template <class T>
class image 
{
public:
    // create an image
    image(const int width, const int height, const bool init = true)
    {
        w = width;
        h = height;
        data = new T[w * h];  // allocate space for image data
        access = new T*[h];   // allocate space for row pointers
        
        // initialize row pointers
        for (int i = 0; i < h; i++)
            access[i] = data + (i * w);  
        
        if (init)
            memset(data, 0, w * h * sizeof(T));
    }

    // delete an image
    ~image()
    {
        delete [] data; 
        delete [] access;
    }

    // copy an image
    image<T> *copy() const
    {      
        image<T> *im = new image<T>(w, h, false);
        memcpy(im->data, data, w * h * sizeof(T));
        return im;
    }
    
    // get the width of an image
    int width() const { return w; }
    
    // get the height of an image
    int height() const { return h; }

    // access image data
    inline T get(int x, int y) const { return access[y][x]; }

    // set image data
    inline void set(int x, int y, T v) { access[y][x] = v; }

    // check if (x,y) is inside the image
    bool check(int x, int y) const { return 0<=x && x<w && 0<=y && y<h; }
  
    T *data;    // image data
    T **access; // row pointers

private:
    int w, h;
};

//-------------------------------------------------------------------------
// color type
//-------------------------------------------------------------------------

typedef struct { unsigned char r, g, b; } rgb;

//-------------------------------------------------------------------------
// image IO
//-------------------------------------------------------------------------

extern void stbi_image_free(void *retval_from_stbi_load);
extern unsigned char *stbi_load(char const *filename, 
								int *x, int *y, int *comp, int req_comp);
static image<rgb>* load(const char* filename, const int bg_color = 255)
{
	int w,h,n;
	unsigned char *data = stbi_load(filename, &w, &h, &n, 0);
	if (data && (n == 3 || n == 4)) {
		image<rgb> *img = new image<rgb>(w, h);
		unsigned char *data_cur = data;
        switch (n) {
        case 3:
            for (int y = 0; y < h; y++) {
                rgb *row_dst = img->access[y]; 
                for (int x = 0; x < w; x++) {
                    row_dst[x].r = *data_cur++;
                    row_dst[x].g = *data_cur++;
                    row_dst[x].b = *data_cur++;
                }
            }
            break;
        case 4:
            for (int y = 0; y < h; y++) {
                rgb *row_dst = img->access[y]; 
                for (int x = 0; x < w; x++) {
                    int r = *data_cur++;
                    int g = *data_cur++;
                    int b = *data_cur++;
                    int a = *data_cur++;
                    int bg = bg_color * (255 - a);
                    row_dst[x].r = (r*a + bg) / 255;
                    row_dst[x].g = (g*a + bg) / 255;
                    row_dst[x].b = (b*a + bg) / 255;
                }
            }
            break;
        default:
            break;
        }
		stbi_image_free(data);
		return img;
	}
	return 0;
}

static void save(const char* filename, const image<rgb>* im)
{
    int width = im->width();
    int height = im->height();
    std::ofstream file(filename, std::ios::out | std::ios::binary);
    file << "P6\n" << width << " " << height << "\n" << UCHAR_MAX << "\n";
    file.write((char *)im->data, width * height * sizeof(rgb));
}