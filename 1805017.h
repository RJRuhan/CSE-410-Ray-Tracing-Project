
#include<bits/stdc++.h>
#include "bitmap_image.hpp"

#ifdef __linux__
#include <GL/glut.h>
#elif WIN32
#include <windows.h>
#include <glut.h>
#endif

using namespace std;

class Object;
class PointLight;
class SpotLight;


#define pi (2*acos(0.0))


extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;
extern vector <Object*> objects;
extern int recursionLevel;

double ephsilon = 1e-6;

class Vector{
public:
    double x, y, z;

    Vector(double x = 0, double y = 0, double z = 0){
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector operator+(Vector v){
        Vector temp;
        temp.x = x + v.x;
        temp.y = y + v.y;
        temp.z = z + v.z;
        return temp;
    }

    Vector operator-(Vector v){
        Vector temp;
        temp.x = x - v.x;
        temp.y = y - v.y;
        temp.z = z - v.z;
        return temp;
    }

    Vector operator-(double a){
        Vector temp;
        temp.x = x - a;
        temp.y = y - a;
        temp.z = z - a;
        return temp;
    }

    Vector operator*(double a){
        Vector temp;
        temp.x = x * a;
        temp.y = y * a;
        temp.z = z * a;
        return temp;
    }

    Vector operator-(){
        Vector temp;
        temp.x = -x;
        temp.y = -y;
        temp.z = -z;
        return temp;
    }

    double operator*(Vector v){
        double a = x*v.x + y*v.y + z*v.z;
        return a;
    }

    double length(){
        return sqrt(x*x + y*y + z*z);
    }


};

extern Vector pos, u, r, l;


Vector crossProduct(Vector a, Vector b) {
    Vector cross;
    cross.x = a.y * b.z - a.z * b.y;
    cross.y = a.z * b.x - a.x * b.z;
    cross.z = a.x * b.y - a.y * b.x;
    return cross;
}


Vector normalize(Vector a) {
    double length = sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    Vector normalized;
    normalized.x = a.x / length;
    normalized.y = a.y / length;
    normalized.z = a.z / length;
    return normalized;
}


Vector rotateAroundPerpAxis(Vector v, Vector axis, double theta) {

    Vector perp = crossProduct(v, axis);

    Vector rotated = v*cos(theta) + perp*sin(theta);
    
    rotated = normalize(rotated);
    return rotated;
}

void moveAlongVector(double step, Vector l) {
    pos.x += l.x * step;
    pos.y += l.y * step;
    pos.z += l.z * step;
}


class PointLight{
public: 
    Vector position;
    double color[3];

    void draw(){
        glPointSize(5);
        glBegin(GL_POINTS);
        glColor3f(color[0], color[1], color[2]);
        glVertex3f(position.x, position.y, position.z);
        glEnd();
    }
};

class SpotLight {
public:
    Vector position;
    double color[3];
    Vector direction;
    double cutOffAngle;

    void draw(){    
        glPointSize(15);
        glBegin(GL_POINTS);
        glColor3f(color[0], color[1], color[2]);
        glVertex3f(position.x, position.y, position.z);
        glEnd();
    }

};

class Ray{
    public:
    Vector start;
    Vector dir;

    Ray(Vector start, Vector dir){
        this->start = start;
        dir = normalize(dir);
        this->dir = dir;

    }

    void draw(){

    }
};



class Object{
public:
    Vector reference_point;
    // should have x, y, z
    double height, width, length;
    double color[3];
    vector<double> coEfficients  = {0,0,0,0}; // default values
    int shine ;  // exponent term of specular component

    Object(){}

    virtual void draw() = 0;

    virtual double intersect(Ray ray, double* color, int level){
        return -1;
    } 

};

class Sphere: public Object{
    public:

    Sphere(){}

    virtual void draw(){

        // cout << "Drawing Sphere" << endl;
        // cout << "Reference Point: " << reference_point.x << " " << reference_point.y << " " << reference_point.z << endl;
        // cout << "Radius: " << length << endl;
        // cout << "Color: " << color[0] << " " << color[1] << " " << color[2] << endl;

        int stacks = 30;
        int slices = 20;

        Vector points[100][100];
        int i, j;
        double h, r;

        for (i = 0; i <= stacks; i++)
        {
            h = length * sin(((double)i / (double)stacks) * (pi / 2));
            r = length * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= slices; j++)
            {
                points[i][j].x = r * cos(((double)j / (double)slices) * 2 * pi);
                points[i][j].y = r * sin(((double)j / (double)slices) * 2 * pi);
                points[i][j].z = h;
            }
        }

        for (i = 0; i < stacks; i++)
        {
            glPushMatrix();
            glTranslatef(reference_point.x, reference_point.y, reference_point.z);
            glColor3f(color[0], color[1], color[2]);
            for (j = 0; j < slices; j++)
            {
                glBegin(GL_QUADS);
                {
                    
                    glVertex3f(points[i][j].x, points[i][j].y, points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, points[i + 1][j].z);
                    
                    glVertex3f(points[i][j].x, points[i][j].y, -points[i][j].z);
                    glVertex3f(points[i][j + 1].x, points[i][j + 1].y, -points[i][j + 1].z);
                    glVertex3f(points[i + 1][j + 1].x, points[i + 1][j + 1].y, -points[i + 1][j + 1].z);
                    glVertex3f(points[i + 1][j].x, points[i + 1][j].y, -points[i + 1][j].z);
                }
                glEnd();
            }
            glPopMatrix();
        }
    }

    double intersectUtility(Ray ray){
        ray.start = ray.start - reference_point; // adjust ray origin
            
        double a = 1;
        double b = 2 * (ray.dir*ray.start);
        double c = (ray.start*ray.start) - (length*length);

        double discriminant = pow(b, 2) - 4 * a * c;
        double t = -1;
        if (discriminant < 0){
            t = -1;
        }
        else{
    
            double t1 = (-b - sqrt(discriminant)) / (2 * a);
            double t2 = (-b + sqrt(discriminant)) / (2 * a);

            if(t1 < 0 && t2 < 0){
                t = -1;
            }
            else if(t1 < 0){
                t = t2;
            }
            else if(t2 < 0){
                t = t1;
            }
            else{
                t = min(t1, t2);
            }
        }

        return t;
    }



    virtual double intersect(Ray ray, double* color, int level){

        double t = intersectUtility(ray);

        if (t < 0){
            return -1;
        }

        if(level == 0){
            return t;
        }

        Vector intersectionPoint = ray.start + ray.dir*t;

        color[0] = this->color[0] * coEfficients[0];
        color[1] = this->color[1] * coEfficients[0];
        color[2] = this->color[2] * coEfficients[0];

        for (int i = 0;i < pointLights.size();i++){
            Vector lightPos = pointLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){

                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*pointLights[i].color[0]*lambert*this->color[0];
                color[1] += coEfficients[1]*pointLights[i].color[1]*lambert*this->color[1];
                color[2] += coEfficients[1]*pointLights[i].color[2]*lambert*this->color[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*pointLights[i].color[0]*specular*this->color[0];
                color[1] += coEfficients[2]*pointLights[i].color[1]*specular*this->color[1];
                color[2] += coEfficients[2]*pointLights[i].color[2]*specular*this->color[2];
            }
        }

        for(int i = 0; i < spotLights.size(); i++){

            Vector lightPos = spotLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;
            lightDir = normalize(lightDir);

            double angle = acos(lightDir*spotLights[i].direction/(lightDir.length()*spotLights[i].direction.length()));
            angle = angle*180/pi;

            if (fabs(angle) >= spotLights[i].cutOffAngle){
                continue;
            }

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){
                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*spotLights[i].color[0]*lambert*this->color[0];
                color[1] += coEfficients[1]*spotLights[i].color[1]*lambert*this->color[1];
                color[2] += coEfficients[1]*spotLights[i].color[2]*lambert*this->color[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*spotLights[i].color[0]*specular*this->color[0];
                color[1] += coEfficients[2]*spotLights[i].color[1]*specular*this->color[1];
                color[2] += coEfficients[2]*spotLights[i].color[2]*specular*this->color[2];
            }
        }

        if(level < recursionLevel){
            Ray normal = getSurfaceNormal(intersectionPoint);
            Ray reflectionRay = Ray(intersectionPoint, ray.dir - normal.dir*2*(ray.dir*normal.dir));
            reflectionRay.start = reflectionRay.start + reflectionRay.dir*ephsilon;
            
            int nearest = -1;
            double tMin = -1,t = -1;
            for(int k = 0; k < objects.size(); k++){
                t = objects[k]->intersect(reflectionRay,color,0);
                if (t > 0 && (nearest == -1 || t < tMin)){
                    tMin = t;
                    nearest = k;
                }
            }

            if(nearest != -1){
                double reflectedColor[3] = {0,0,0};
                t = objects[nearest]->intersect(reflectionRay, reflectedColor, level+1);
                
                color[0] += reflectedColor[0]*coEfficients[3];
                color[1] += reflectedColor[1]*coEfficients[3];
                color[2] += reflectedColor[2]*coEfficients[3];
                
            }
        }

        return t;
    }

    Ray getSurfaceNormal(Vector point){        
        return Ray(point, point - reference_point);
    }
};

double determinant(double mat[3][3]){
    double d1 = mat[0][0] * (mat[1][1]*mat[2][2] - mat[1][2]*mat[2][1]);
    double d2 = mat[0][1] * (mat[1][0]*mat[2][2] - mat[1][2]*mat[2][0]);
    double d3 = mat[0][2] * (mat[1][0]*mat[2][1] - mat[1][1]*mat[2][0]);
    return d1 - d2 + d3;
}

class Triangle: public Object{
    public:
    Vector a, b, c;
    Triangle(){}

    virtual void draw(){
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(a.x, a.y, a.z);
            glVertex3f(b.x, b.y, b.z);
            glVertex3f(c.x, c.y, c.z);
        }
        glEnd();
    }

    double intersectUtility(Ray ray){

        double AMat[3][3] {
            {a.x - b.x, a.x - c.x, ray.dir.x},
            {a.y - b.y, a.y - c.y, ray.dir.y},
            {a.z - b.z, a.z - c.z, ray.dir.z}
        };

       double betaMat[3][3] = {
            {a.x - ray.start.x, a.x - c.x, ray.dir.x},
            {a.y - ray.start.y, a.y - c.y, ray.dir.y},
            {a.z - ray.start.z, a.z - c.z, ray.dir.z}
        };

        double gammaMat[3][3] = {
            {a.x - b.x, a.x - ray.start.x, ray.dir.x},
            {a.y - b.y, a.y - ray.start.y, ray.dir.y},
            {a.z - b.z, a.z - ray.start.z, ray.dir.z}
        };

        double tMat[3][3] = {
            {a.x - b.x, a.x - c.x, a.x - ray.start.x},
            {a.y - b.y, a.y - c.y, a.y - ray.start.y},
            {a.z - b.z, a.z - c.z, a.z - ray.start.z}
        };
        

        double Adet = determinant(AMat);
        double beta = determinant(betaMat) / Adet;
        double gamma = determinant(gammaMat) / Adet;
        double t = determinant(tMat) / Adet;

        if (beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0){
            return t;
        }
        else{
            return -1;
        }
    }

    virtual double intersect(Ray ray, double* color, int level){

        double t = intersectUtility(ray);

        if (t < 0){
            return -1;
        }

        if(level == 0){
            return t;
        }

        Vector intersectionPoint = ray.start + ray.dir*t;

        color[0] = this->color[0] * coEfficients[0];
        color[1] = this->color[1] * coEfficients[0];
        color[2] = this->color[2] * coEfficients[0];

        for (int i = 0;i < pointLights.size();i++){
            Vector lightPos = pointLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint, lightRay);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){

                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*pointLights[i].color[0]*lambert*this->color[0];
                color[1] += coEfficients[1]*pointLights[i].color[1]*lambert*this->color[1];
                color[2] += coEfficients[1]*pointLights[i].color[2]*lambert*this->color[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*pointLights[i].color[0]*specular*this->color[0];
                color[1] += coEfficients[2]*pointLights[i].color[1]*specular*this->color[1];
                color[2] += coEfficients[2]*pointLights[i].color[2]*specular*this->color[2];
            }
        }

        for(int i = 0; i < spotLights.size(); i++){

            Vector lightPos = spotLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;
            lightDir = normalize(lightDir);

            double angle = acos(lightDir*spotLights[i].direction/(lightDir.length()*spotLights[i].direction.length()));
            angle = angle*180/pi;

            if (fabs(angle) >= spotLights[i].cutOffAngle){
                continue;
            }

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint,lightRay);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){
                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*spotLights[i].color[0]*lambert*this->color[0];
                color[1] += coEfficients[1]*spotLights[i].color[1]*lambert*this->color[1];
                color[2] += coEfficients[1]*spotLights[i].color[2]*lambert*this->color[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*spotLights[i].color[0]*specular*this->color[0];
                color[1] += coEfficients[2]*spotLights[i].color[1]*specular*this->color[1];
                color[2] += coEfficients[2]*spotLights[i].color[2]*specular*this->color[2];
            }
        }

        if(level < recursionLevel){

            Ray normal = getSurfaceNormal(intersectionPoint,ray);
            Ray reflectionRay = Ray(intersectionPoint, ray.dir - normal.dir*2*(ray.dir*normal.dir));
            reflectionRay.start = reflectionRay.start + reflectionRay.dir*ephsilon;
            
            int nearest = -1;
            double tMin = -1,t = -1;
            for(int k = 0; k < objects.size(); k++){
                t = objects[k]->intersect(reflectionRay,color,0);
                if (t > 0 && (nearest == -1 || t < tMin)){
                    tMin = t;
                    nearest = k;
                }
            }

            if(nearest != -1){
                double reflectedColor[3] = {0,0,0};
                t = objects[nearest]->intersect(reflectionRay, reflectedColor, level+1);
                
                color[0] += reflectedColor[0]*coEfficients[3];
                color[1] += reflectedColor[1]*coEfficients[3];
                color[2] += reflectedColor[2]*coEfficients[3];
                
            }
        }

        return t;
    }

    Ray getSurfaceNormal(Vector point,Ray incidentRay){
        Vector normal = crossProduct(b - a, c - a);
        normal = normalize(normal);

        if (normal*incidentRay.dir < 0){
            normal = -normal;
        }

        return Ray(point, normal);
        
    }
};

class GeneralQuadraticSurface: public Object{
    public:
    double A, B, C, D, E, F, G, H, I, J;

    GeneralQuadraticSurface(){}

    virtual void draw(){

    }

    bool validPoint(Vector point){
        if ( fabs(length) > ephsilon && ( point.x < reference_point.x || point.x > reference_point.x + length ) ) return false;
        if ( fabs(width) > ephsilon && ( point.y < reference_point.y || point.y > reference_point.y + width ) ) return false;
        if ( fabs(height) > ephsilon && ( point.z < reference_point.z || point.z > reference_point.z + height ) ) return false;
        return true;
    }

    virtual double intersectUtility(Ray ray){

        double X0 = ray.start.x;
        double Y0 = ray.start.y;
        double Z0 = ray.start.z;

        double X1 = ray.dir.x;
        double Y1 = ray.dir.y;
        double Z1 = ray.dir.z;

        double C0 = A*X1*X1 + B*Y1*Y1 + C*Z1*Z1 + D*X1*Y1 + E*X1*Z1 + F*Y1*Z1;
        double C1 = 2*A*X0*X1 + 2*B*Y0*Y1 + 2*C*Z0*Z1 + D*(X0*Y1 + X1*Y0) + E*(X0*Z1 + X1*Z0) + F*(Y0*Z1 + Y1*Z0) + G*X1 + H*Y1 + I*Z1;
        double C2 = A*X0*X0 + B*Y0*Y0 + C*Z0*Z0 + D*X0*Y0 + E*X0*Z0 + F*Y0*Z0 + G*X0 + H*Y0 + I*Z0 + J;

        double discriminant = C1*C1 - 4*C0*C2;

        if(discriminant < 0) return -1;

        if(fabs(C0) < 1e-5) {
            return -C2/C1;
        }
        double t1 = (-C1 - sqrt(discriminant))/(2*C0);
        double t2 = (-C1 + sqrt(discriminant))/(2*C0);

        double t = -1;

        if(t1 < 0 && t2 < 0){
            t = -1;
        }
        else if(t1 < 0){
            t = t2;
            Vector intersectionPoint = ray.start + ray.dir*t;
            if(validPoint(intersectionPoint)){
                return t;
            }
        }
        else if(t2 < 0){
            t = t1;
            Vector intersectionPoint = ray.start + ray.dir*t;
            if(validPoint(intersectionPoint)){
                return t;
            }
        }
        else{
            if(t2<t1) swap(t1,t2);
            t = t1;
            Vector intersectionPoint = ray.start + ray.dir*t;
            if(validPoint(intersectionPoint)){
                return t;
            }

            t = t2;
            intersectionPoint = ray.start + ray.dir*t;
            if(validPoint(intersectionPoint)){
                return t;
            }
        }

        return -1;
    }

    virtual double intersect(Ray ray, double* color, int level){

        double t = intersectUtility(ray);

        if (t < 0){
            return -1;
        }

        if(level == 0){
            return t;
        }

        Vector intersectionPoint = ray.start + ray.dir*t;

        color[0] = this->color[0] * coEfficients[0];
        color[1] = this->color[1] * coEfficients[0];
        color[2] = this->color[2] * coEfficients[0];

        for (int i = 0;i < pointLights.size();i++){
            Vector lightPos = pointLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;
            lightDir = normalize(lightDir);

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){

                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*pointLights[i].color[0]*lambert*this->color[0];
                color[1] += coEfficients[1]*pointLights[i].color[1]*lambert*this->color[1];
                color[2] += coEfficients[1]*pointLights[i].color[2]*lambert*this->color[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*pointLights[i].color[0]*specular*this->color[0];
                color[1] += coEfficients[2]*pointLights[i].color[1]*specular*this->color[1];
                color[2] += coEfficients[2]*pointLights[i].color[2]*specular*this->color[2];
            }
        }

        for(int i = 0; i < spotLights.size(); i++){

            Vector lightPos = spotLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;

            double angle = acos(lightDir*spotLights[i].direction/(lightDir.length()*spotLights[i].direction.length()));
            angle = angle*180/pi;

            if (fabs(angle) >= spotLights[i].cutOffAngle){
                continue;
            }

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){
                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*spotLights[i].color[0]*lambert*this->color[0];
                color[1] += coEfficients[1]*spotLights[i].color[1]*lambert*this->color[1];
                color[2] += coEfficients[1]*spotLights[i].color[2]*lambert*this->color[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*spotLights[i].color[0]*specular*this->color[0];
                color[1] += coEfficients[2]*spotLights[i].color[1]*specular*this->color[1];
                color[2] += coEfficients[2]*spotLights[i].color[2]*specular*this->color[2];
            }
        }

        if(level < recursionLevel){
            Ray normal = getSurfaceNormal(intersectionPoint);
            Ray reflectionRay = Ray(intersectionPoint, ray.dir - normal.dir*2*(ray.dir*normal.dir));
            reflectionRay.start = reflectionRay.start + reflectionRay.dir*ephsilon;
            
            int nearest = -1;
            double tMin = -1,t = -1;
            for(int k = 0; k < objects.size(); k++){
                t = objects[k]->intersect(reflectionRay,color,0);
                if (t > 0 && (nearest == -1 || t < tMin)){
                    tMin = t;
                    nearest = k;
                }
            }

            if(nearest != -1){
                double reflectedColor[3] = {0,0,0};
                t = objects[nearest]->intersect(reflectionRay, reflectedColor, level+1);
                
                color[0] += reflectedColor[0]*coEfficients[3];
                color[1] += reflectedColor[1]*coEfficients[3];
                color[2] += reflectedColor[2]*coEfficients[3];
                
            }
        }

        return t;
    }



    Ray getSurfaceNormal(Vector point){
        double X = point.x;
        double Y = point.y;
        double Z = point.z;

        double X0 = 2*A*X + D*Y + E*Z + G;
        double Y0 = 2*B*Y + D*X + F*Z + H;
        double Z0 = 2*C*Z + E*X + F*Y + I;

        Vector normal = {X0, Y0, Z0};
        return Ray(point, normal);
        
    }

};

class Floor: public Object{
    
    public:
    int tileCount;


    Floor(double FloorWidth, double TileWidth){
        reference_point.x = -FloorWidth/2;
        reference_point.y = -FloorWidth/2;
        reference_point.z = 0;
        width = FloorWidth;
        length = TileWidth;
        height = 0;
        tileCount = FloorWidth/TileWidth;
    }

    virtual void draw(){
        for (int i = 0; i < tileCount; i++)
		{
			for (int j = 0; j < tileCount; j++)
			{
				if (((i + j) % 2) == 0) 
                    glColor3f(1, 1, 1);
				else 
                    glColor3f(0, 0, 0);

				glBegin(GL_QUADS);
				{
					glVertex3f(reference_point.x + i * length, reference_point.y + j * length, 0);
					glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + j * length, 0);
					glVertex3f(reference_point.x + (i + 1) * length, reference_point.y + (j + 1) * length, 0);
					glVertex3f(reference_point.x + i * length, reference_point.y + (j + 1) * length, 0);
				}
				glEnd();
			}
		}
    }

    void getColorAtPoint(Vector point,double* color){

        int tileX = (int)(point.x-reference_point.x)/length;
        int tileY = (int)(point.y-reference_point.y)/length;

        if( tileX < 0 || tileX >= tileCount || tileY < 0 || tileY >= tileCount){
            color[0] = 1;
            color[1] = 1;
            color[2] = 1;
            return;
        }

        if((tileX+tileY)%2 == 0){
            color[0] = 1;
            color[1] = 1;
            color[2] = 1;
        }
        else{
            color[0] = 0;
            color[1] = 0;
            color[2] = 0;
        }

    }

    double intersectUtility(Ray ray){
        if (fabs(ray.dir.z) < ephsilon){
            return -1;
        }
        double t = -ray.start.z/ray.dir.z;

        Vector intersectionPoint = ray.start + ray.dir*t;
        if (intersectionPoint.x <= reference_point.x || intersectionPoint.x > -reference_point.x || intersectionPoint.y <= reference_point.y || intersectionPoint.y > -reference_point.y){
            return -1;
        }

        return t;
    }

    virtual double intersect(Ray ray, double* color, int level){

        double t = intersectUtility(ray);

        if (t < 0){
            return -1;
        }

        if(level == 0){
            return t;
        }

        Vector intersectionPoint = ray.start + ray.dir*t;

        double tileColor[3] = {0,0,0};
        getColorAtPoint(intersectionPoint,tileColor);

        color[0] = tileColor[0] * coEfficients[0];
        color[1] = tileColor[1] * coEfficients[0];
        color[2] = tileColor[2] * coEfficients[0];

        for (int i = 0;i < pointLights.size();i++){
            Vector lightPos = pointLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){

                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*pointLights[i].color[0]*lambert*tileColor[0];
                color[1] += coEfficients[1]*pointLights[i].color[1]*lambert*tileColor[1];
                color[2] += coEfficients[1]*pointLights[i].color[2]*lambert*tileColor[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*pointLights[i].color[0]*specular*tileColor[0];
                color[1] += coEfficients[2]*pointLights[i].color[1]*specular*tileColor[1];
                color[2] += coEfficients[2]*pointLights[i].color[2]*specular*tileColor[2];
            }
        }

        for(int i = 0; i < spotLights.size(); i++){

            Vector lightPos = spotLights[i].position;
            Vector lightDir = intersectionPoint - lightPos;
            lightDir = normalize(lightDir);

            double angle = acos(lightDir*spotLights[i].direction/(lightDir.length()*spotLights[i].direction.length()));
            angle = angle*180/pi;

            if (fabs(angle) >= spotLights[i].cutOffAngle){
                continue;
            }

            Ray lightRay = Ray(lightPos, lightDir);

            Ray normal = getSurfaceNormal(intersectionPoint);

            Vector distanceVector = intersectionPoint - lightPos;
            double distance = distanceVector.length();
            if(distance < ephsilon) 
                continue;

            bool rayBlocked = false;
            for (int j = 0;j < objects.size();j++){
                double t = objects[j]->intersect(lightRay, color, 0);
                if(t > 0 && distance - t > ephsilon){
                    rayBlocked = true;
                    break;
                }
            }

            if(!rayBlocked){
                double lambert = max(0.0, -lightRay.dir*normal.dir);
                color[0] += coEfficients[1]*spotLights[i].color[0]*lambert*tileColor[0];
                color[1] += coEfficients[1]*spotLights[i].color[1]*lambert*tileColor[1];
                color[2] += coEfficients[1]*spotLights[i].color[2]*lambert*tileColor[2];

                Ray reflection = Ray(intersectionPoint, lightRay.dir - normal.dir*2*(lightRay.dir*normal.dir));
                double specular = pow(max(0.0, -ray.dir*reflection.dir), shine);
                
                color[0] += coEfficients[2]*spotLights[i].color[0]*specular*tileColor[0];
                color[1] += coEfficients[2]*spotLights[i].color[1]*specular*tileColor[1];
                color[2] += coEfficients[2]*spotLights[i].color[2]*specular*tileColor[2];
            }
        }

        if(level < recursionLevel){
            Ray normal = getSurfaceNormal(intersectionPoint);
            Ray reflectionRay = Ray(intersectionPoint, ray.dir - normal.dir*2*(ray.dir*normal.dir));
            reflectionRay.start = reflectionRay.start + reflectionRay.dir*ephsilon;
            
            int nearest = -1;
            double tMin = -1,t = -1;
            for(int k = 0; k < objects.size(); k++){
                t = objects[k]->intersect(reflectionRay,color,0);
                if (t > 0 && (nearest == -1 || t < tMin)){
                    tMin = t;
                    nearest = k;
                }
            }

            if(nearest != -1){
                double reflectedColor[3] = {0,0,0};
                t = objects[nearest]->intersect(reflectionRay, reflectedColor, level+1);
                
                color[0] += reflectedColor[0]*coEfficients[3];
                color[1] += reflectedColor[1]*coEfficients[3];
                color[2] += reflectedColor[2]*coEfficients[3];
                
            }
        }

        return t;
    }

    Ray getSurfaceNormal(Vector point){
            return Ray(point, Vector(0, 0, 1));
    }
};
