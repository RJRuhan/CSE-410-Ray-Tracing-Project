#include <bits/stdc++.h>
using namespace std;
#include "1805017.h"



Vector pos, u, r, l;

vector<Object*> objects;
vector<PointLight> pointLights;
vector<SpotLight> spotLights;

int recursionLevel;
int imageWidth, imageHeight;

void display()
{
    glEnable(GL_DEPTH_TEST);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    

    gluLookAt(
        pos.x, pos.y, pos.z,
        pos.x + l.x, pos.y + l.y, pos.z + l.z,
        u.x, u.y, u.z
        );

    cout << "pos: " << pos.x << " " << pos.y << " " << pos.z << endl;
    cout << "l: " << l.x << " " << l.y << " " << l.z << endl;
    cout << "u: " << u.x << " " << u.y << " " << u.z << endl;
    cout << "r: " << r.x << " " << r.y << " " << r.z << endl;

    for (int i = 0; i < objects.size(); i++){
        objects[i]->draw();
    }

    for (int i = 0; i < pointLights.size(); i++){
        pointLights[i].draw();
    }

    for (int i = 0; i < spotLights.size(); i++){
        spotLights[i].draw();
    }

    glutSwapBuffers();
}

int imageCount = 1;
double windowWidth = 500, windowHeight = 500;
double viewAngle = 80;

void capture(){
    bitmap_image image = bitmap_image(imageWidth, imageHeight);
    image.set_all_channels(0, 0, 0);

    // for(int i=0;i<imageWidth;i++)
	// 	for(int j=0;j<imageHeight;j++)
	// 		image.set_pixel(i, j, 0, 0, 0);

    double planeDistance = (windowHeight / 2.0) / tan((viewAngle / 2.0) * (pi / 180.0));
    Vector topLeft = pos + l * planeDistance - r * (windowWidth / 2.0) + u * (windowHeight / 2.0);

    double du = windowWidth / (double)imageWidth;
    double dv = windowHeight / (double)imageHeight;

    topLeft = topLeft + r * (du * 0.5) - u * (dv * 0.5);

    int nearest = -1;
    double t,tMin;

    for (int i = 0;i < imageWidth;i++){
        for(int j = 0;j < imageHeight;j++){
            Vector currPixel = topLeft + r * du * i - u * dv * j;

            Ray ray(pos, currPixel - pos);
            double color[3] = {0,0,0};

            tMin = -1;
            nearest = -1;
            for(int k = 0; k < objects.size(); k++){
                t = objects[k]->intersect(ray,color,0);
                if (t > 0 && (nearest == -1 || t < tMin)){
                    tMin = t;
                    nearest = k;
                }
            }

            if (nearest != -1){
                
                color[0] = 0, color[1] = 0, color[2] = 0;

                double t = objects[nearest]->intersect(ray,color,1);

                if(color[0] > 1){
                    color[0] = 1;
                    cout << "Color Red Exceeded 1" << endl;
                }

                if(color[1] > 1){
                    color[1] = 1;
                    cout << "Color Green Exceeded 1" << endl;
                }

                if(color[2] > 1){
                    color[2] = 1;
                    cout << "Color Blue Exceeded 1" << endl;
                }

                if(color[0] < 0){
                    color[0] = 0;
                    cout << "Color Red Less than 0" << endl;
                }

                if(color[1] < 0){
                    color[1] = 0;
                    cout << "Color Green Less than 0" << endl;
                }

                if(color[2] < 0){
                    color[2] = 0;
                    cout << "Color Blue Less than 0" << endl;
                }

                image.set_pixel(i, j, color[0] * 255, color[1] * 255, color[2] * 255);

            }
        }
    }

    image.save_image("Output_1"+to_string(imageCount)+".bmp");
	imageCount++;
    cout << "Image Captured" << endl;
}



void keyboardHandler(unsigned char key, int x, int y){
    cout << key << endl;
    switch(key){
        case '0':
			capture();
			break;
        case '1':
            cout << "1 pressed" << endl;
            l = rotateAroundPerpAxis(l,u,-0.1);
            r = rotateAroundPerpAxis(r,u,-0.1);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case '2':
            cout << "2 pressed" << endl;
            l = rotateAroundPerpAxis(l,u,0.1);
            r = rotateAroundPerpAxis(r,u,0.1);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case '3':
            cout << "3 pressed" << endl;
            l = rotateAroundPerpAxis(l,r,0.1);
            u = rotateAroundPerpAxis(u,r,0.1);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case '4':
            cout << "4 pressed" << endl;
            l = rotateAroundPerpAxis(l,r,-0.1);
            u = rotateAroundPerpAxis(u,r,-0.1);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case '5':
            cout << "5 pressed" << endl;
            u = rotateAroundPerpAxis(u,l,0.1);
            r = rotateAroundPerpAxis(r,l,0.1);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case '6':
            cout << "6 pressed" << endl;
            u = rotateAroundPerpAxis(u,l,-0.1);
            r = rotateAroundPerpAxis(r,l,-0.1);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        
    }
}


void specialKeyHandler(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
            cout << "Up arrow key pressed" << endl;
            moveAlongVector(1.0,l);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case GLUT_KEY_DOWN:
            cout << "Down arrow key pressed" << endl;
            moveAlongVector(-1.0,l);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case GLUT_KEY_LEFT:
            cout << "Left arrow key pressed" << endl;
            moveAlongVector(-1.0,r);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case GLUT_KEY_RIGHT:
            cout << "Right arrow key pressed" << endl;
            moveAlongVector(1.0,r);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case GLUT_KEY_PAGE_UP:
            cout << "Page up key pressed" << endl;
            moveAlongVector(1.0,u);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
        case GLUT_KEY_PAGE_DOWN:
            cout << "Page down key pressed" << endl;
            moveAlongVector(-1.0,u);
            glutPostRedisplay(); // Trigger display function to update the scene
            break;
            
    }

}

void loadData(){

    ifstream input("scene.txt");

    input >> recursionLevel;
    input >> imageWidth;
    imageHeight = imageWidth;

    int objectCount;
    input >> objectCount;

    for (int i = 0;i < objectCount;i++){
        string type;
        input >> type;

        if (type == "sphere"){
            Sphere* s = new Sphere();
            input >> s->reference_point.x >> s->reference_point.y >> s->reference_point.z;
            input >> s->length;
            input >> s->color[0] >> s->color[1] >> s->color[2];
            input >> s->coEfficients[0] >> s->coEfficients[1] >> s->coEfficients[2] >> s->coEfficients[3];
            input >> s->shine;
            objects.push_back(s);
        }
        else if (type == "triangle"){
            Triangle* t = new Triangle();
            input >> t->a.x >> t->a.y >> t->a.z;
            input >> t->b.x >> t->b.y >> t->b.z;
            input >> t->c.x >> t->c.y >> t->c.z;
            input >> t->color[0] >> t->color[1] >> t->color[2];
            input >> t->coEfficients[0] >> t->coEfficients[1] >> t->coEfficients[2] >> t->coEfficients[3];
            input >> t->shine;
            objects.push_back(t);
        }
        else if (type == "general"){
            GeneralQuadraticSurface* g = new GeneralQuadraticSurface();
            input >> g->A >> g->B >> g->C >> g->D >> g->E >> g->F >> g->G >> g->H >> g->I >> g->J;
            input >> g->reference_point.x >> g->reference_point.y >> g->reference_point.z;
            input >> g->width >> g->length >> g->height;
            input >> g->color[0] >> g->color[1] >> g->color[2];
            input >> g->coEfficients[0] >> g->coEfficients[1] >> g->coEfficients[2] >> g->coEfficients[3];
            input >> g->shine;
            objects.push_back(g);
        }
    }

    int pointLightCount;
    input >> pointLightCount;

    for(int i = 0;i < pointLightCount;i++){
        PointLight p;
        input >> p.position.x >> p.position.y >> p.position.z;
        input >> p.color[0] >> p.color[1] >> p.color[2];
        pointLights.push_back(p);
    }

    int spotLightCount;
    input >> spotLightCount;

    for(int i = 0;i < spotLightCount;i++){
        SpotLight s;
        input >> s.position.x >> s.position.y >> s.position.z;
        input >> s.color[0] >> s.color[1] >> s.color[2];
        input >> s.direction.x >> s.direction.y >> s.direction.z;
        input >> s.cutOffAngle;
        spotLights.push_back(s);
    }

    input.close();
}


void createFloor(int floorWidth,int tileWidth){
    Floor* f = new Floor(floorWidth, tileWidth);
    f->coEfficients = {0.4, 0.2, 0.2, 0.2};
    f->color[0] = 0.5;
    f->color[1] = 0.5;
    f->color[2] = 0.5;
    objects.push_back(f);
}


// pos: 102.119 9.73766 17.1955
// l: -0.999364 -0.0317246 -0.0162548
// u: -0.0220075 0.190379 0.981464
// r: -0.028042 0.981198 -0.190956

// pos: -44.6606 -90.4065 35.1252
// l: 0.476714 0.877026 -0.0597398
// u: -0.051503 0.0957079 0.994076
// r: 0.877548 -0.470813 0.090794

// 200,0,10
// 0,0,1
// -1 / sqrt(2), 1 / sqrt(2), 0
// -1 / sqrt(2), -1 / sqrt(2), 0
void init()
{   
    pos = Vector(5, 5, 5);
    u = Vector(0,0,1);
    r = Vector(-0.707107, 0.707107, 0);
    l = Vector(-0.707107, -0.707107, 0);

    loadData();

    createFloor(1000, 20);

    glMatrixMode(GL_PROJECTION); 
    glLoadIdentity();
    gluPerspective(viewAngle, 1, 1, 1000);

}

template <typename T>
void deleteAndClearVector(std::vector<T*>& vec) {
    for (auto& element : vec) {
        delete element;  // Delete the pointer
    }

    vec.clear();  // Clear the vector (shrink to size 0)
    vec.shrink_to_fit();  // Shrink the vector to fit
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Offline 2: Ray Tracing");

    init();

    glutDisplayFunc(display);

    glutKeyboardFunc(keyboardHandler);
    glutSpecialFunc(specialKeyHandler);
    glutMainLoop();

    objects.erase(objects.begin(), objects.end());

    deleteAndClearVector(objects);

    pointLights.clear();
    pointLights.shrink_to_fit();
    
    spotLights.clear();
    spotLights.shrink_to_fit();

    return 0;
}