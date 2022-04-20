 #version 330 

uniform vec4 eye;
uniform vec4 ambient;
uniform vec4[20] objects;
uniform vec4[20] objColors;
uniform vec4[10] lightsDirection;
uniform vec4[10] lightsIntensity;
uniform vec4[10] lightsPosition;
uniform ivec4 sizes;

in vec3 position0;
in vec3 normal0;

layout (location = 0) out vec4 outColor;

float intersection(inout int sourceIndx,vec3 sourcePoint,vec3 v)
{
    float tmin = 1.0e10;
    int indx = -1;
    for(int i=0;i<sizes.x;i++) //every object
    {
        if(i==sourceIndx)
            continue;
        if(objects[i].w > 0) //sphere
        {
            vec3 p0o =  objects[i].xyz - sourcePoint;
            float r = objects[i].w;
            float b = dot(v,p0o);
            float delta = b*b - dot(p0o,p0o) + r*r;
             float t;
            if(delta >= 0)
            {
                if(b>=0)
                    t = b - sqrt(delta);
                else
                    t = b + sqrt(delta);
                if(t<tmin && t>0)
                {
                    tmin = t;
                    indx = i;
                } 
            }   
        }
        else  //plane
        {    
            vec3 n =  normalize(objects[i].xyz);
            vec3 p0o = -objects[i].w*n/length(objects[i].xyz) - sourcePoint;
            float t = dot(n,p0o)/dot(n,v); 
            if(t>0 && t<tmin)
            {
                tmin = t;
                indx = i;
            }
        }
    }
    sourceIndx = indx; 
    return tmin;
}


//body index in objects, point on surface of object, diffuseFactor for plane squares
vec3 colorCalc(int sourceIndx, vec3 sourcePoint,vec3 u,float diffuseFactor)
{
    vec3 color = ambient.rgb*objColors[sourceIndx].rgb;
    float specularCoeff = 0.7f;
    for(int i = 0;i<sizes.y;i++) //every light source
    {
        vec3 v;
        if(lightsDirection[i].w < 0.5 ) //directional
        {
            int indx = sourceIndx;
            v = normalize(lightsDirection[i].xyz);
           //  v = normalize(vec3(0.0,0.5,-1.0));
            float t = intersection(indx,sourcePoint,-v);
            
            // TODO: tamir, why??? planes are see through?
            if(indx < 0 || objects[indx].w<=0) //no intersection
             {
               // vec3 u = normalize(sourcePoint - eye.xyz);
                if(objects[sourceIndx].w > 0) //sphere
                {
                    
                    vec3 n = -normalize( sourcePoint - objects[sourceIndx].xyz);
                    vec3 refl = normalize(reflect(v,n));
                    if(dot(v,n)>0.0 )
                        color+= max(specularCoeff * lightsIntensity[i].rgb * pow(dot(refl,u),objColors[sourceIndx].a),vec3(0.0,0.0,0.0));  //specular  
                    color+= max(diffuseFactor * objColors[sourceIndx].rgb * lightsIntensity[i].rgb * dot(v,n),vec3(0.0,0.0,0.0)) ;  //difuse
                    //        color = vec3(1.0,1.0,0.0);
                }
                else  //plane
                {
                    vec3 n = normalize(objects[sourceIndx].xyz);
                    vec3 refl = normalize(reflect(v,n));
                    
                    color = min(color + max(specularCoeff * lightsIntensity[i].rgb * pow(dot(refl,u),objColors[sourceIndx].a),vec3(0.0,0.0,0.0)),vec3(1.0,1.0,1.0)); //specular
                    color = min( color + max(diffuseFactor * objColors[sourceIndx].rgb * lightsIntensity[i].rgb * dot(v,n),vec3(0.0,0.0,0.0)),vec3(1.0,1.0,1.0)); //difuse
                 
                  //      color = vec3(1.0,1.0,0.0);
                }
            }
         //   else if(indx == 1)
          //          color = lightsIntensity[i].rgb;
            
        }
        else  //flashlight
        {
            int indx = -1;
            v = -normalize(lightsPosition[i].xyz - sourcePoint);
            if(dot(v,normalize(lightsDirection[i].xyz))<lightsPosition[i].w)
            {
                continue;
            }
            else
            {
                //vec3 u = normalize(sourcePoint - eye.xyz);
                float t = intersection(indx,lightsPosition[i].xyz,v);
                if(indx == sourceIndx) //no intersection
                {
                    if(objects[sourceIndx].w > 0) //sphere
                    {
                        vec3 n = -normalize( sourcePoint - objects[sourceIndx].xyz);
                        vec3 refl = normalize(reflect(v,n));
                        if(dot(v,n)>0.0)
                          color+=max(specularCoeff * lightsIntensity[i].rgb * pow(dot(refl,u),objColors[sourceIndx].a),vec3(0.0,0.0,0.0)); //specular
                        color+= max(diffuseFactor * objColors[sourceIndx].rgb * lightsIntensity[i].rgb * dot(v,n),vec3(0.0,0.0,0.0));
                      //          color = vec3(1.0,1.0,0.0);            
                    }
                    else  //plane
                    {

                        vec3 n = normalize(objects[sourceIndx].xyz);
                        vec3 refl = normalize(reflect(v,n)); //specular
                        color = min(color + max(specularCoeff * lightsIntensity[i].rgb * pow(dot(refl,u),objColors[sourceIndx].a),vec3(0.0,0.0,0.0)),vec3(1.0,1.0,1.0));
                        color = min(color + max(diffuseFactor * objColors[sourceIndx].rgb * lightsIntensity[i].rgb *dot(v,n),vec3(0.0,0.0,0.0)),vec3(1.0,1.0,1.0));
                       // color = vec3(1.0,1.0,0.0);
                    }
                }
                //else if(indx == 1)
                //    color = lightsIntensity[i].rgb;
            }
        }
    }
         //   color = vec3(1.0,1.0,0.0);
    return min(color,vec3(1.0,1.0,1.0));
}

void findIntersection(out float dist, out vec3 normal, out vec3 intersectionPoint, vec4 object, vec3 p0, vec3 ray) {
    dist = -1.0;
    normal = vec3(-1);
    intersectionPoint = vec3(-1);
    if(object.w <= 0) {
        //plane

        float d = object.w;
        normal = normalize(object.xyz);
        dist = -(dot(normal, p0) + d) / dot(normal, ray);
        intersectionPoint = p0 + ray*dist;
    }
    else {
        //sphere
        vec3 o = object.xyz;
        float r = object.w;
        vec3 L = o - p0;
        float tm = dot(L, ray);
        float dSquared = pow(length(L), 2) - pow(tm, 2);
        if(dSquared <= pow(r, 2)) {
            float th = sqrt(pow(r, 2) - dSquared);
            float t1 = tm - th, t2 = tm + th;
            if(t1 > 0) {
                if(t2 > 0) {
                    dist = min(t1, t2);
                }
                else {
                    dist = t1;
                }
            }
            else {
                dist = t2;
            }
            intersectionPoint = p0 + ray*dist;
            normal = normalize(intersectionPoint - o);
        }
    }
}

void findFirstIntersectingObject(out int intersectionIndex, out float intersectionDistance, out vec3 intersectionPoint, out vec3 intersectionNormal, vec3 p0, vec3 ray) {
    float minDist = -1;
    int minObjectIndex = -1;
    vec3 minInterPoint = vec3(-1);
    vec3 minInterNormal = vec3(-1);
    for(int i = 0; i < sizes[0]; i++) {
        vec4 curObject = objects[i];
        float dist;
        vec3 interNormal, interPoint;
        findIntersection(dist, interNormal, interPoint, curObject, p0, ray);
        if(dist > 1.5e-6 && (dist < minDist || minObjectIndex == -1)) {
            minDist = dist;
            minObjectIndex = i;
            minInterPoint = interPoint;
            minInterNormal = interNormal;
        }
    }
    intersectionIndex = minObjectIndex;
    intersectionDistance = minDist;
    intersectionPoint = minInterPoint;
    intersectionNormal = minInterNormal;
}

float squareMod(float c, float factor, float sideLength) {
    return mod(int((factor * c) / sideLength), 2);
}

struct PlaneMap {
    vec3 b1;
    vec3 b2;
    vec3 anchor;
};

PlaneMap mapPlaneToXyPlane(vec4 plane) {
    // See https://stackoverflow.com/questions/8780646/mapping-coordinates-from-plane-given-by-normal-vector-to-xy-plane.

    vec3 planeNormal = plane.xyz;
    float d = plane.w;

    vec3 b1, b2;
    if (plane.z == 0) {
        if (plane.y == 0) {
            float x = -d / plane.x;
            b1 = vec3(x, 0, 1);
            b2 = vec3(x, 1, 0);
        }
        else if (plane.x == 0) {
            float y = -d / plane.y;
            b1 = vec3(1, y, 0);
            b2 = vec3(0, y, 1);
        }
        else {
            // a*x + b*y + d = 0 -->
            // y = -(a/b)*x - d/b

            // get the perpendicular line in a 2D plane which goes through the origin
            float m = plane.y / plane.x;
            // find the intersection point between the above line and `y = m*x`
            float x = -(d/plane.y) / (m + 1/m);
            float y = m * x;
            b1 = vec3(x, y, 1);
            vec2 b2_2d = vec2(x, y) + normalize(vec2(-plane.y, plane.x));
            b2 = vec3(b2_2d, 0);
        }
    }
    else {
        b1 = vec3(1, 0, -plane.y / plane.z);
        b2 = vec3(0, 1, -plane.x / plane.z);
    }
    b1 = normalize(b1);
    b2 = normalize(b2);

    vec3 anchor = (d/dot(planeNormal, planeNormal)) * planeNormal;
    return PlaneMap(b1, b2, anchor);
}

vec3 applyPlaneMap(PlaneMap map, vec3 point) {
    mat2x3 coefficients = mat2x3(map.b1, map.b2);
    mat3x2 coefficientsT = transpose(coefficients);
    mat3x2 leftInverse = inverse(coefficientsT * coefficients) * coefficientsT;
    return vec3(leftInverse * point, 0);
}

bool isDarkSquare(vec3 point) {
    return (mod(int(1.5 * point.x), 2) == mod(int(1.5 * point.y), 2)) == ((point.x < 0) == (point.y < 0));
}

vec3 calculateColor_noTracing(vec3 vRay, vec3 point, vec3 pointNormal, int objectIndex) {
    vec3 objectColor = objColors[objectIndex].xyz;
    float shinniness = objColors[objectIndex].w;
    vec4 object = objects[objectIndex];
    vec3 color = objectColor * ambient.xyz;
    
    vec3 specularFactors = vec3(0.7);
    vec3 diffuseFactors = vec3(1);
    if (object.w <= 0) {
        float squareSideLength = 1 / 1.5;
        float d = object.w;
        diffuseFactors = vec3(1);
        PlaneMap map = mapPlaneToXyPlane(object);
        if (isDarkSquare(applyPlaneMap(map, point))) {
            diffuseFactors = vec3(0.5);
        }
    }

    int currentSpotlightIdx = -1;
    for(int i = 0; i < sizes[1]; i++) {
        vec4 curLight = lightsDirection[i];
        vec3 lightDirection = curLight.xyz;
        vec4 lightIntensity = lightsIntensity[i];

        vec3 vPointToLight;
        vec3 vPointToLightUnnormalized;
        float cosBetween;
        vec3 intensity = lightIntensity.xyz;
//        intensity *= dot(vPointToLight, pointNormal);
        if(curLight.w < 0.5) {
            //directional
            vPointToLightUnnormalized = -lightDirection;
            vPointToLight = normalize(vPointToLightUnnormalized);
            lightDirection = normalize(lightDirection);
            // TODO: why without this it looks bad?
            intensity *= dot(vPointToLight, pointNormal);
        }
        else {
            //spotlight
            currentSpotlightIdx += 1;
            vec4 spotlightInfo = lightsPosition[currentSpotlightIdx];
            vec3 spotlightPosition = spotlightInfo.xyz;
            float spotlightHalfApertureCos = spotlightInfo.w;
            lightDirection = normalize(lightDirection);
            vPointToLightUnnormalized = spotlightPosition - point;
            vPointToLight = normalize(vPointToLightUnnormalized);
            float cosBetween = dot(-vPointToLight, lightDirection);
            if(cosBetween <= spotlightHalfApertureCos) {
                // in range
                // TODO: do we need to?
                intensity *= cosBetween;
            }
            else {
                continue;
            }
        }

        int blockingObject;
        float blockingDist;
        vec3 blockingPoint, blockingNormal;
        findFirstIntersectingObject(blockingObject, blockingDist, blockingPoint, blockingNormal, point, vPointToLight);
        if (
            blockingObject >= 0
//            && blockingObject != objectIndex
            && ((blockingObject == objectIndex) != (objects[blockingObject].w > 0)) // A plane cannot block itself and a sphere can
            && (curLight.w < 0.5 || blockingDist < length(vPointToLightUnnormalized)) // Obly block when either if the light is directional or if the 'blocking' object appears between the point and the spotlight
//            && objects[blockingObject].w > 0
        ) {
            continue;
        }

        // TODO: why this max?
        vec3 diffuse = diffuseFactors * objectColor * intensity * dot(pointNormal, vPointToLight);
        diffuse = max(diffuse, vec3(0));
        vec3 refl = normalize(reflect(-vPointToLight, pointNormal));
        vec3 specular = specularFactors * intensity * pow(dot(-vRay, refl), shinniness);
        // TODO: why this max?
        specular = max(specular, vec3(0));
        if(object.w > 0) {
            //sphere

            // TODO: why this if?
            //       and why is the `vPointToLight` should not be negated?
            if (dot(vPointToLight, pointNormal) > 0) {
                color += specular;
            }
            color += diffuse;
        }
        else {
            // plane

            // TODO: why this min?
            color = min(color + specular, vec3(1));
            color = min(color + diffuse, vec3(1));
        }
    }

    return color;
}

void main()
{  
    vec3 vRay = normalize(position0.xyz - eye.xyz);
    int interObject;
    float interDist;
    vec3 interPoint, interNormal;
    findFirstIntersectingObject(interObject, interDist, interPoint, interNormal, position0, vRay);
//    interObject = -1;
//    float t = intersection(interObject, position0, vRay);
//    interPoint = position0 + t*vRay;


    vec3 color;
    if(interObject == -1) {
        color = vec3(1, 1, 1);
    }
    else {
        color = calculateColor_noTracing(vRay, interPoint, interNormal, interObject);
//        float diffuseFactor = 1;
//        if (objects[interObject].w <= 0 && isDarkSquare(interPoint)) {
//            diffuseFactor = 0.5;
//        }
//        color = colorCalc(interObject, interPoint, vRay, diffuseFactor), 1;
//        color = vec4(1, 0, 0, 1);
    }
    outColor = vec4(color, 1);
}
 


