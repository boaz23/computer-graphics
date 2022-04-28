#version 330 

uniform vec4 eye;
uniform vec4 ambient;
uniform vec4[20] objects;
uniform vec4[20] objColors;
uniform vec4[10] lightsDirection;
uniform vec4[10] lightsIntensity;
uniform vec4[10] lightsPosition;
uniform ivec4 sizes; //{number of objects , number of lights , number of reflective objects, 0}  
uniform vec4 screenCenter;

in vec3 position0;


#define delta 0.00000000001

#define inf 1.0/0.0
#define vecinf vec4(inf, inf, inf, -1.0)

bool isPlane(vec4 object)
{
    return object.w < 0.0;
}

bool isSpotlight(vec4 lightDir)
{
    return lightDir.w > 0.0;
}

vec3 normal(vec4 object, vec3 P){
    return !isPlane(object) ? normalize(P - object.xyz) :  -normalize(object.xyz);
}

vec3 getDirectionVector()
{
    vec3 P = position0;
    vec3 tv = P - eye.xyz;
    return normalize(tv);
}

vec3 getDirectionVector(vec3 p0, vec3 p1)
{
    vec3 v = normalize(p1 - p0);
    return v;
}


float planeIntersection(vec3 P0, vec3 v, vec4 plane)
{
    vec3 N = plane.xyz;
    float a = plane.x;
    float b = plane.y;
    float c = plane.z;
    float d = plane.w;
    
    //calculating a P on the plane
    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    if (b != 0.0)
        y = -d/b;
    else if (a != 0.0)
        x = -d/a;
    else if (c != 0.0)
        z = -d/c;
    vec3 Q0 = vec3(x,y,z);

    vec3 over = Q0 - P0;
    float under = dot(N, v);
    if(abs(under) < delta)
    {
        return inf;
    }

    float t = dot(N, over)/under;
    if (t < 0.0)
    {
        return inf;
    }
    return t;
}

float sphereIntersection(vec3 P0, vec3 v, vec4 sphere)
{
    vec3 O = sphere.xyz;
    vec3 L = O - P0;
    float tm = dot(L,v);
    float d2 = dot(L,L) - tm*tm;
    float r = sphere.w;
    if(d2 > r*r)
    {
        return inf;
    }
    float th = sqrt(r*r - d2);
    float t;
    if ((tm - th) > 0.0)
    {
        t = tm - th;
    }
    else if ((tm + th) > 0.0)
    {
        t = tm + th;
    }
    else
    {
        t = inf;
    }
    return t;
}

float getIntersection(vec3 P0, vec3 v, vec4 object)
{
    if(isPlane(object))
    {
        return planeIntersection(P0, v, object);
    }
    else
    {
        return sphereIntersection(P0, v, object);
    }
    return inf;
}

//vec2 intersection(vec3 P0,vec3 v, float targetDist, int thisNdx, bool shadow)
//{
//    float closestObject = -1.0;
//    float closestT = 0.0;
//    float closestDist = targetDist;
//    float t;
//    vec3 P;
//    float dist;
//    for(int i = 0; i < sizes[0]; ++i){
//        if(i != thisNdx && !(shadow && isPlane(objects[i]))){
//            t = getIntersection(P0, v, objects[i]);
//            P = P0 + t*v;
//            if(t < inf  && (dist = distance(P0, P)) < closestDist)
//            {   
//                closestDist = dist;
//                closestT = t;
//                closestObject = float(i);
//            }
//        }
//    }
//    return vec2(closestT, closestObject);
//}

struct Intersection {
    int objectIndex;
    float distance;
    vec3 point;
    vec3 pointNormal;
};

Intersection findIntersection(vec4 object, vec3 p0, vec3 ray) {
    float dist = -1.0;
    vec3 normal = vec3(-1);
    vec3 intersectionPoint = vec3(-1);
    if (isPlane(object)) {
        //plane
        float d = object.w;
        vec3 perpendicular = object.xyz;
        normal = normalize(perpendicular);
        dist = -(dot(perpendicular, p0) + d) / dot(perpendicular, ray);
        intersectionPoint = p0 + ray*dist;
    }
    else {
        //sphere
        vec3 o = object.xyz;
        float r = object.w;
        vec3 L = o - p0;
        float tm = dot(L, ray);
        float dSquared = pow(length(L), 2) - pow(tm, 2);
        if (dSquared <= pow(r, 2)) {
            float th = sqrt(pow(r, 2) - dSquared);
            float t1 = tm - th, t2 = tm + th;
            if (t1 >= delta) {
                dist = t1;
            }
            else if (t2 >= delta) {
                dist = t2;
            }
            else {
                dist = -1.0;
            }
            intersectionPoint = p0 + ray*dist;
            normal = normalize(intersectionPoint - o);
        }
    }
    return Intersection(-1, dist, intersectionPoint, normal);
}

vec2 intersection(vec3 P0,vec3 v, float targetDist, int thisNdx, bool shadow) {
    Intersection minIntersection = Intersection(-1, -1, vec3(-1), vec3(-1));
    for(int i = 0; i < sizes[0]; i++) {
        vec4 curObject = objects[i];

        // Skip transparent planes. The light is considered to go through them without change
        if (false
            || i == thisNdx
            || (shadow && isPlane(curObject))
        ) {
            continue;
        }

        Intersection intersection = findIntersection(curObject, P0, v);
        float t = intersection.distance;
        if (true
            && t >= delta && t < targetDist && t < inf
            && (minIntersection.objectIndex == -1 || t < minIntersection.distance)
        ) {
            intersection.objectIndex = i;
            minIntersection = intersection;
            if (shadow) {
                break;
            }
        }
    }

    return vec2(minIntersection.distance, float(minIntersection.objectIndex));
}

vec3 defuse(vec3 N, vec3 L, vec3 Kd){
    float cosTh = dot(normalize(N), normalize(L));
    vec3 Idef = Kd * cosTh;
    return clamp(Idef, 0.0, 1.0);
}

float ourPow(float x, float e) {
    int e_i = int(e);
    float p = 1.0;
    for (int i = 0; i < e; i++) {
        p *= x;
    }
    return p;
}

vec3 specular(vec3 negV,vec3 N, vec3 L, float n){
    // from wikipedia https://en.wikipedia.org/wiki/Phong_reflection_model
//    vec3 R = 2*(dot(normalize(L)  , N)) * N - normalize(L);
    vec3 R = -reflect(L, N);
    float RDotV = max(0, dot(negV , R));
    return clamp(vec3(0.7,0.7,0.7) *  pow(RDotV, n),0.0, 1.0);
}

vec3 getGridColor(vec3 P, vec4 objColor){
    if(mod(int(1.5*P.x + 150), 2) == mod(int(1.5*P.y + 150), 2))
        return vec3(0.5 , 0.5 , 0.5) * objColor.xyz;
    else
        return objColor.xyz;
}

vec3 defAndSpec(vec3 P, vec4 lightPos, vec4 lightInt, vec4 lightDir, vec4 object, int objNdx,  vec4 objColor, vec3 v)
{
    vec3 N = normal(object, P);
    vec3 Kd  = objColor.xyz;
    if(isPlane(object)){
        Kd = getGridColor(P, objColor);
    }

    float n = objColor.w;

    vec3 L;

    bool shadow = false;
    vec2 blocking = vec2(inf, 1.0);
    if(!isSpotlight(lightDir)){
        blocking = intersection(P, L, inf, objNdx, true);
        L = -normalize(lightDir.xyz);
    }
    else{
        float cosdeg = dot(normalize(lightDir.xyz), getDirectionVector(lightPos.xyz, P));
        if(acos(cosdeg) < (acos(lightPos.w))){
            L = getDirectionVector(P, lightPos.xyz);
            blocking = intersection(lightPos.xyz, -L, inf, objNdx, true);
        }
        else {
            blocking = vec2(inf, 1.0);
        }
    }
    
    shadow = (int(blocking.y) != -1);
    
    if (!shadow){
        vec3 def = defuse(N, L,  Kd);
//        vec3 def = vec3(0);
        vec3 spec = specular(-v,  N , L, n);
        return (def+spec)*lightInt.xyz;
    }
    return vec3(0.0, 0.0, 0.0);
}

vec3 colorCalc(vec3 P0, float t, vec4 object, int objNdx, vec4 objColor , vec3 v)
{
    vec3 P = P0 + t*v;
    vec3 Ka = objColor.xyz; 
    vec3 Ia = ambient.xyz;
    vec3 I = vec3(0.0, 0.0, 0.0);
    I += clamp(Ka*Ia, 0.0, 1.0);
    int spotNdx = 0;
    for(int i = 0; i <sizes[1]; ++i)
    {   
        if(!isSpotlight(lightsDirection[i]))
        {
            I += defAndSpec(P, vecinf, lightsIntensity[i], lightsDirection[i], object, objNdx, objColor, v);
        }
        else
        {
            I += defAndSpec(P, lightsPosition[spotNdx], lightsIntensity[i], lightsDirection[i], object, objNdx, objColor, v);
            spotNdx++;
        }
    }
    return I;
}

vec3 reflect_vec(vec3 u , vec3 n){
    vec3 u_n = -n * (dot(u ,n));
    vec3 u_star = u + u_n;
    vec3 v = u_n + u_star; //  v = u_n + u^* = u + 2u_n = u -2n(dot(u,n))
    return normalize(v);
}

vec3 reflection(vec3 P0, float t, vec4 object, int objNdx, vec4 objColor, vec3 v)
{
    vec3 P = P0;
    vec3 N;
    vec2 t_and_objNdx;
    int firstObjNdx = objNdx;
    float firstT = t;
    vec3 firstV = v;
    for(int i = 0; i <=4; ++i){
        P = P + t*v;
        N = normal(object, P);
        v = reflect(v, N);
        t_and_objNdx = intersection(P, v, inf, objNdx, false);
        t = t_and_objNdx.x;
        objNdx = int(t_and_objNdx.y);
        object = objects[objNdx];
        objColor = objColors[objNdx];
        if(objNdx >= int(sizes[2]))
        {
            break;
        }
        else if(objNdx == -1)
        {
            return vec3(0.0,0.0,0.0);
        }
    }
    return colorCalc(P, t, object, objNdx, objColor, v);
}

void main()
{    
    vec3 v = getDirectionVector();
    vec2 t_and_objNdx = intersection(eye.xyz, v, inf, -1, false);
    float t = t_and_objNdx.x;
    int objNdx = int(t_and_objNdx.y);
    vec3 color = vec3(0.0, 1.0, 0.0);
    if(t > 0.0)
    {
        if(float(objNdx) >= sizes[2])
            color = colorCalc(eye.xyz, t, objects[objNdx], objNdx, objColors[objNdx], v);
        else
            color = reflection(eye.xyz, t, objects[objNdx], objNdx, objColors[objNdx], v);
    }
    gl_FragColor = vec4(clamp(color, 0.0, 1.0), 1.0);
}