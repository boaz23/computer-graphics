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

void findIntersection(out float dist, out vec3 normal, out vec3 intersectionPoint, vec4 object, vec3 p0, vec3 ray){
    dist = -1.0;
    normal = vec3(-1);
    intersectionPoint = vec3(-1);
    if(object.w <= 0){
        //plane
        normal = normalize(object.xyz);
        dist = -(dot(normal, p0) + object.w) / dot(normal, ray);
        intersectionPoint = p0 + ray*dist;
    }
    else{
        //sphere
        vec3 L = object.xyz - p0;
        float tm = dot(L, ray);
        float dSquared = pow(length(L), 2) - pow(tm, 2);
        if(dSquared <= pow(object.w, 2)){
            float th = sqrt(pow(object.w, 2) - dSquared);
            if(tm - th > 0){
                if(tm + th > 0){
                    dist = min(tm - th, tm + th);
                }
                else{
                    dist = tm - th;
                }
            }
            else{
                dist = tm + th;
            }
            intersectionPoint = p0 + ray*dist;
            normal = normalize(intersectionPoint - object.xyz);
        }
    }
}

void findFirstIntersectingObject(out int intersectionIndex, out float intersectionDistance, out vec3 intersectionPoint, out vec3 intersectionNormal, vec3 p0, vec3 ray){
    float minDist = -1;
    int minObjectIndex = -1;
    vec3 minInterPoint = vec3(-1);
    vec3 minInterNormal = vec3(-1);
    for(int i = 0; i<sizes[0]; i++){
        vec4 curObject = objects[i];
        float dist;
        vec3 interNormal, interPoint;
        findIntersection(dist, interNormal, interPoint, curObject, p0, ray);
        if(dist > 0.01 && (dist < minDist || minObjectIndex == -1)){
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

void main()
{  
    vec3 vRay = normalize(position0.xyz - eye.xyz);
    int interObject;
    float interDist;
    vec3 interPoint, interNormal;
    findFirstIntersectingObject(interObject, interDist, interPoint, interNormal, eye.xyz, vRay);
    if(interObject == -1){
        gl_FragColor = vec4(1, 0, 0, 1);
    }
    else{
        vec4 color = objColors[interObject] * ambient;
        int currentSpotlightIdx = 0;
        for(int i = 0; i < sizes[1]; i++){
            vec4 curLight = lightsDirection[i];
            vec3 rayToLight;
            float cosBetween;
            vec4 intensity;
            if(curLight.w < 0.5){
                //directional
                rayToLight = -curLight.xyz;
                intensity = lightsIntensity[i] * dot(normalize(rayToLight), interNormal);
            }
            else{
                //spotlight
                rayToLight = lightsPosition[currentSpotlightIdx].xyz - interPoint;
                float cosBetween = dot(normalize(-rayToLight), normalize(curLight.xyz));
                if(cosBetween <= lightsPosition[currentSpotlightIdx].w){
                    // in range
                    intensity = lightsIntensity[i] * cosBetween;
                }
                else{
                    currentSpotlightIdx += 1;
                    continue;
                }
                currentSpotlightIdx += 1;
            }
            int blockingObject;
            float blockingDist;
            vec3 blockingPoint, blockingNormal;
            findFirstIntersectingObject(blockingObject, blockingDist, blockingPoint, blockingNormal, interPoint, rayToLight);
            if(blockingDist > 0 && blockingObject != interObject && (curLight.w < 0.5 || blockingDist < length(rayToLight))){
                continue;
            }
            float x = interPoint.x;
            float y = interPoint.y;
/*            if(objects[interObject].w <= 0 && (((mod(int(1.5*x),2) == mod(int(1.5*y),2)) && ((x>0 && y>0) || (x<0 && y<0))) || ((mod(int(1.5*x),2) != mod(int(1.5*y),2) && ((x<0 && y>0) || (x>0 && y<0))))))
                gl_FragColor = vec4(colorCalc(indx,p,v,0.5),1);
        else 
            gl_FragColor = vec4(colorCalc(indx,p,v,1.0),1);      
            */
            color += objColors[interObject] * dot(interNormal, normalize(rayToLight)) * intensity;
        }
        gl_FragColor = color;
    }
}
 


