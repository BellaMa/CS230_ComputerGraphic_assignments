/**
 * ray_tracer.cpp
 * CS230
 * -------------------------------
 * Implement ray tracer here.
 */

#define SET_RED(P, C)   (P = (((P) & 0x00ffffff) | ((C) << 24)))
#define SET_GREEN(P, C)  (P = (((P) & 0xff00ffff) | ((C) << 16)))
#define SET_BLUE(P, C) (P = (((P) & 0xffff00ff) | ((C) << 8)))

#include "ray_tracer.h"

using namespace std;

const double Object::small_t=1e-6;
//--------------------------------------------------------------------------------
// utility functions
//--------------------------------------------------------------------------------
double sqr(const double x)
{
    return x*x;
}
//typedef unsigned int Pixel;
Pixel Pixel_Color(const Vector_3D<double>& color)
{
    Pixel pixel=0;
    SET_RED(pixel,(unsigned char)(min(color.x,1.0)*255));
    SET_GREEN(pixel,(unsigned char)(min(color.y,1.0)*255));
    SET_BLUE(pixel,(unsigned char)(min(color.z,1.0)*255));
    return pixel;
}
//--------------------------------------------------------------------------------
// Shader
//--------------------------------------------------------------------------------
Vector_3D<double> Phong_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color, light_dir,tem, emitted_light,h,v,r;
    Vector_3D<double> I_ambiend, I_diffuse, I_specular;
    // TODO: determine the color

    v = tem - ray.direction;

    
    for(int i = 0; i < world.lights.size(); ++i)
    {
        // ambient
        emitted_light = world.lights.at(i)->Emitted_Light(ray);
        I_ambiend = color_ambient * emitted_light;

        light_dir = world.lights.at(i)->position - intersection_point;
        light_dir.Normalize();

        // add a small value to the starting point to avoid dark spots
        Ray shadow(intersection_point + same_side_normal * intersection_object.small_t, light_dir);

        
        // diffuse
        double tem_l_n;
        tem_l_n = tem.Dot_Product(light_dir,same_side_normal);
        I_diffuse =  color_diffuse * emitted_light * max(tem_l_n,0.0);
        
        // specular
        double tem_v_r;
        r = same_side_normal * 2.0 * tem.Dot_Product(same_side_normal,light_dir)  - light_dir;
        r.Normalize();
        tem_v_r = tem.Dot_Product(v,r);
        I_specular = color_specular * emitted_light * pow(max(0.0, tem_v_r),specular_power);

     
        if(world.Closest_Intersection(shadow) != 0 && world.enable_shadows)
            color += I_ambiend ;
        else
            color += I_ambiend + I_diffuse + I_specular;
     
   } 
            
  

    return color;
}

Vector_3D<double> Reflective_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    Vector_3D<double> color, light_dir,tem, emitted_light,h,v,r;
    Vector_3D<double> I_ambiend, I_diffuse, I_specular;
    // TODO: determine the color

    v = tem - ray.direction;

    
    for(int i = 0; i < world.lights.size(); ++i)
    {
        // ambient
        emitted_light = world.lights.at(i)->Emitted_Light(ray);
        I_ambiend = color_ambient * emitted_light;

        light_dir = world.lights.at(i)->position - intersection_point;
        light_dir.Normalize();

        Ray shadow(intersection_point+ same_side_normal * intersection_object.small_t, light_dir);
        world.Closest_Intersection(shadow);


        
        // diffuse
        double tem_l_n;
        tem_l_n = tem.Dot_Product(light_dir,same_side_normal);
        I_diffuse =  color_diffuse * emitted_light * max(tem_l_n,0.0);
        
        // specular
        double tem_v_r;
        r = same_side_normal * 2.0 * tem.Dot_Product(same_side_normal,light_dir)  - light_dir;
        r.Normalize();
        tem_v_r = tem.Dot_Product(v,r);
        I_specular = color_specular * emitted_light * pow(max(0.0, tem_v_r),specular_power);


        
        if(world.Closest_Intersection(shadow) != 0  && world.enable_shadows)
            color += I_ambiend ;
        else
            color += I_ambiend + I_diffuse + I_specular;
     
    }

    double tem_n_v = tem.Dot_Product(same_side_normal,v);
    Vector_3D<double> reflec = same_side_normal*2.0*tem_n_v - v;

    Ray reflection(intersection_point + same_side_normal*intersection_object.small_t,reflec);

    reflection.recursion_depth = ray.recursion_depth + 1;
    color += world.Cast_Ray(reflection,ray) * reflectivity;



    return color;
}

Vector_3D<double> Flat_Shader::
Shade_Surface(const Ray& ray,const Object& intersection_object,const Vector_3D<double>& intersection_point,const Vector_3D<double>& same_side_normal) const
{
    return color;
}

//--------------------------------------------------------------------------------
// Objects
//--------------------------------------------------------------------------------
// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Sphere::
Intersection(Ray& ray) const
{
    // TODO
    double delta,bb,ac;
    Vector_3D<double> tem;
    Vector_3D<double> e_c = ray.endpoint - center;
    bb = 2.0 * tem.Dot_Product(ray.direction,e_c);
    ac = e_c.Length_Squared() - sqr(radius);
    delta = sqr(bb) - 4.0 * ac;


    if(delta < 0)
        return false; // no intersection

    double t1 = (-bb - sqrt(delta)) / 2;
    if(t1 > 0)
    {
        if(ray.t_max == 0)
        {
            // ray hasn't hit anything
            ray.semi_infinite = false;
            ray.t_max = t1;
            ray.current_object = this;
            return true;
        }
        if(t1 < ray.t_max)
        {
            // although ray has hit something, but this obj is closer
            ray.t_max = t1;
            ray.current_object = this;
            return true;
        }
    }

    // there's something closer or this obj is at opposite dirction
    return true;
}

Vector_3D<double> Sphere::
Normal(const Vector_3D<double>& location) const
{
    Vector_3D<double> normal;
    // TODO: set the normal
    normal = location - center;
    normal.Normalize();
    return normal;
}

// determine if the ray intersects with the sphere
// if there is an intersection, set t_max, current_object, and semi_infinite as appropriate and return true
bool Plane::
Intersection(Ray& ray) const
{
    // TODO
    Vector_3D<double> tem_x1_x = x1 - ray.endpoint;
    double den, num, tem_t;
    Vector_3D<double> tem;
    den = tem.Dot_Product(tem_x1_x,normal);
    num = tem.Dot_Product(ray.direction,normal);
    tem_t = den / num;
    if(tem_t > 0)
    {
        if(!ray.semi_infinite)
        {
            if(ray.t_max > tem_t)
            {
                ray.t_max = tem_t;
                ray.current_object = this;
                
                return true;
            }
        }else
        {
            ray.t_max = tem_t;
            ray.current_object = this;
            ray.semi_infinite = false;

            return true;
        }
    }

    return false;
}

Vector_3D<double> Plane::
Normal(const Vector_3D<double>& location) const
{
    return normal;
}
//--------------------------------------------------------------------------------
// Camera
//--------------------------------------------------------------------------------
// Find the world position of the input pixel
Vector_3D<double> Camera::
World_Position(const Vector_2D<int>& pixel_index)
{
    Vector_3D<double> result;
    // TODO
    result = focal_point + horizontal_vector*(film.pixel_grid.X(pixel_index).x);
    result += vertical_vector*(film.pixel_grid.X(pixel_index).y);
    return result;
}
//--------------------------------------------------------------------------------
// Render_World
//--------------------------------------------------------------------------------
// Find the closest object of intersection and return a pointer to it
//   if the ray intersects with an object, then ray.t_max, ray.current_object, and ray.semi_infinite will be set appropriately
//   if there is no intersection do not modify the ray and return 0
const Object* Render_World::
Closest_Intersection(Ray& ray)
{
    // TODO
    for(int i = 0; i < objects.size(); ++i)
        objects.at(i)->Intersection(ray);  

    return ray.current_object;
}

// set up the initial view ray and call
void Render_World::
Render_Pixel(const Vector_2D<int>& pixel_index)
{
    // TODO
    Ray ray; // TODO: set up the initial view ray here
    ray.endpoint = camera.position;
    ray.direction = camera.World_Position(pixel_index) - camera.position;
    ray.direction.Normalize();
    
    Ray dummy_root;
    Vector_3D<double> color=Cast_Ray(ray,dummy_root);
    camera.film.Set_Pixel(pixel_index,Pixel_Color(color));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
Vector_3D<double> Render_World::
Cast_Ray(Ray& ray,const Ray& parent_ray)
{
    // TODO
    Vector_3D<double> color;
    Vector_3D<double> p_intersection, n;
    const Object * curr = this->Closest_Intersection(ray);

    if(curr == 0 ) 
        return color; 
    else if(ray.recursion_depth > this->recursion_depth_limit)
        return color;
    else 
    {
        p_intersection = ray.Point(ray.t_max);
        n = curr->Normal(p_intersection);
        return curr->material_shader->Shade_Surface(ray, *curr, p_intersection, n);
    }

}
