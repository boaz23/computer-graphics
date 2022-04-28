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

#define EPSILON 1e-10

#define clampColor(color) (clamp((color), 0, 1))

#define isZeroP(x) (x < EPSILON)
#define isZero(x) isZeroP(abs(x))

#define OBJ_KIND_PLANE 0
#define OBJ_KIND_SPHERE 1
#define isObjectOfKind(object, kind) (((object).w <= 0) == ((kind) == OBJ_KIND_PLANE))
#define getObjectKind(object) (isObjectOfKind(object, OBJ_KIND_PLANE) ? OBJ_KIND_PLANE : OBJ_KIND_SPHERE)

#define OBJ_FLAGS_TRANSPARENT (1 << 0)
#define OBJ_FLAGS_REFLECTIVE (1 << 1)
int getObjectFlags(int i) {
    if (i < 0) {
        return 0;
    }
    else if (i < sizes[3]) {
        return OBJ_FLAGS_TRANSPARENT;
    }
    else if (i < sizes[3] + sizes[2]) {
        return OBJ_FLAGS_REFLECTIVE;
    }
    else {
        return 0;
    }
}

struct Object {
    int kind;
    vec4 info;
    vec3 color;
    float shinniness;
};

Object getObject(int i) {
    vec4 info = objects[i];
    vec4 colorInfo = objColors[i];

    int kind = getObjectKind(info);
    vec3 color = colorInfo.xyz;
    float shinniness = colorInfo.w;

    return Object(kind, info, color, shinniness);
}

bool isTransparentPlane(int i) {
    vec4 object = objects[i];
    return true
        && isObjectOfKind(object, OBJ_KIND_PLANE)
        && (getObjectFlags(i) & OBJ_FLAGS_TRANSPARENT) != 0;
}

#define LIGHT_KIND_DIRECTIONAL 0
#define LIGHT_KIND_SPOTLIGHT 1
struct Light {
    int kind;
    vec3 intensity;
    vec3 direction;

    // spotlight only
    vec3 position;
    float spotlight_halfApertureCos;
};

Light getLight(int i, int positionIndex) {
    vec4 light = lightsDirection[i];
    vec4 positionInfo;
    int kind;
    if (light.w < 0.5) {
        kind = LIGHT_KIND_DIRECTIONAL;
        positionInfo = vec4(vec3(0), 1);
    }
    else {
        kind = LIGHT_KIND_SPOTLIGHT;
        positionInfo = lightsPosition[positionIndex];
    }

    vec3 direction = normalize(light.xyz);
    vec3 intensity = lightsIntensity[i].xyz;
    vec3 position = positionInfo.xyz;
    float spotlight_halfApertureCos = positionInfo.w;
    return Light(kind, intensity, direction, position, spotlight_halfApertureCos);
}

struct StraightLine {
    vec3 p0;
    vec3 v;
};

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
    if (isObjectOfKind(object, OBJ_KIND_PLANE)) {
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
        float dSquared = (length(L)*length(L)) - (tm*tm);
        if (dSquared <= r*r) {
            float th = sqrt(r*r - dSquared);
            float t1 = tm - th, t2 = tm + th;
            if (!isZeroP(t1)) {
                dist = t1;
            }
            else if (!isZeroP(t2)) {
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

Intersection findFirstIntersectingObject(vec3 p0, vec3 ray, bool skipTransparentPlanes, int skipObjectIndex, bool isForShadowing) {
    Intersection minIntersection = Intersection(-1, -1, vec3(-1), vec3(-1));
    for(int i = 0; i < sizes[0]; i++) {
        vec4 curObject = objects[i];

        // Skip transparent planes. The light is considered to go through them without change
        if (false
            || i == skipObjectIndex
            || (skipTransparentPlanes && isTransparentPlane(i))
            || (isForShadowing &&isObjectOfKind(curObject, OBJ_KIND_PLANE))
        ) {
            continue;
        }

        Intersection intersection = findIntersection(curObject, p0, ray);
        if (true
            && !isZeroP(intersection.distance)
            && (minIntersection.objectIndex == -1 || intersection.distance < minIntersection.distance)
        ) {
            intersection.objectIndex = i;
            minIntersection = intersection;
            if (isForShadowing) {
                break;
            }
        }
    }
    return minIntersection;
}

struct BasedPlane {
    vec3 anchor;
    vec3 base1;
    vec3 base2;
};

BasedPlane findBaseForPlane(vec4 plane) {
    float d = plane.w;
    vec3 planeCoefficients = plane.xyz;
    vec3 planeNormal = normalize(planeCoefficients);
    vec3 anchor = -(d / dot(planeCoefficients, planeCoefficients)) * planeCoefficients;
    vec3 b1, b2;
    if (plane.x == 0 && plane.y == 0) {
        b1 = vec3(1, 0, 0);
        b2 = vec3(0, 1, 0);
    }
    else {
        b1 = normalize(vec3(-plane.y, plane.x, 0));
        b2 = normalize(cross(b1, planeNormal));
    }
    return BasedPlane(anchor, b1, b2);
}

vec2 calculateBaseCoordinates(BasedPlane plane, vec3 point) {
    // `mat2x3` is actually a 3x2 matrix.
    // Solve the following linear equation system:
    //   anchor + mat2x3(base1, base2) * vec2(x, y) = point
    // Equivalently, solve:
    //   mat2x3(base1, base2) * vec2(x, y) = point - anchor
    // The solution vec2(x, y)
    mat2x3 coefficients = mat2x3(plane.base1, plane.base2);
    mat3x2 coefficientsT = transpose(coefficients);
    mat3x2 leftInverse = inverse(coefficientsT * coefficients) * coefficientsT;
    return leftInverse * (point - plane.anchor);
}

#define PLANE_SQUARE_SIDE_LENGTH (1 / 1.5)
int calcSquareIndex(float x) {
    return int(mod(int(x / PLANE_SQUARE_SIDE_LENGTH), 2));
}

bool isDarkSquare(vec2 point) {
    bool equalSqaureIndex = calcSquareIndex(point.x) == calcSquareIndex(point.y);
    // need to specify both `<=` and `>=` because otherwise the colors in the quaters
    // of the plane where their sign differs will be the opposite of what it should be.
    bool equalSign = ((point.x >= 0) == (point.y >= 0)) || ((point.x <= 0) == (point.y <= 0));
    return equalSqaureIndex == equalSign;
}

float calculateDiffuseFactor(vec4 object, vec3 point) {
    float diffuseFactor =1;
    if (isObjectOfKind(object, OBJ_KIND_PLANE)) {
        BasedPlane plane = findBaseForPlane(object);
        vec2 mappedCoordinates = calculateBaseCoordinates(plane, point);
        if (isDarkSquare(mappedCoordinates)) {
            diffuseFactor = 0.5;
        }
    }

    return diffuseFactor;
}

bool isShadowing(Intersection hit, Intersection blocking, Light light) {
    bool isValid = blocking.objectIndex >= 0;
    bool isSameObject = blocking.objectIndex == hit.objectIndex;
    bool isSamePoint = (hit.point != blocking.point && isZeroP(length(hit.point - blocking.point)));
    
    vec4 object = objects[hit.objectIndex];
    bool isObjectASphere = isObjectOfKind(object, OBJ_KIND_SPHERE);

    vec4 blockingObject = objects[blocking.objectIndex];
    bool isBlockingObjectASphere = isObjectOfKind(blockingObject, OBJ_KIND_SPHERE);
    
    return isValid;
    // return isValid && (!isSameObject || false);
    // TODO: why are planes transparent?
    // return isValid && (!isSameObject || false) && isBlockingObjectASphere;
    // TODO: why does this makes weird black points on the spheres?
    // return isValid && (!isSameObject || (isObjectASphere && isSamePoint));
}

vec3 calculateColor_noTracing(vec3 vRay, Intersection intersection) {
    Object object = getObject(intersection.objectIndex);

    vec3 specularFactors = vec3(0.7);
    vec3 diffuseFactors = vec3(calculateDiffuseFactor(object.info, intersection.point));

    vec3 color = object.color;
    color *= ambient.xyz;

    int spotlightIndex = 0;
    for(int i = 0; i < sizes[1]; i++) {
        Light light = getLight(i, spotlightIndex);

        vec3 vPointToLight;
        vec3 intensity = light.intensity;
        Intersection blockingIntersection;
        if (light.kind == LIGHT_KIND_DIRECTIONAL) {
            vPointToLight = -light.direction;

            blockingIntersection = findFirstIntersectingObject(
                intersection.point,
                vPointToLight,
                true,
                intersection.objectIndex,
                true
            );
        }
        else {
            spotlightIndex += 1;
            vPointToLight = normalize(light.position - intersection.point);
            float cosBetween = dot(-vPointToLight, light.direction);
            if (cosBetween < light.spotlight_halfApertureCos) {
                continue;
            }
            else {
                // in range
                
                intensity *= cosBetween;
                // By shooting the ray from the light position in the direction to the point, we avoid
                // hitting objects which are out of the spotlight's range.
                blockingIntersection = findFirstIntersectingObject(
                    light.position,
                    -vPointToLight,
                    true,
                    intersection.objectIndex,
                    true
                );
            }
        }

        if (isShadowing(intersection, blockingIntersection, light)) {
            continue;
        }

        float cosIncoming = dot(intersection.pointNormal, vPointToLight);
        vec3 normal = intersection.pointNormal;
        vec3 diffuse = diffuseFactors * object.color * intensity * cosIncoming;
        vec3 refl = normalize(reflect(-vPointToLight, normal));
        vec3 specular = specularFactors * intensity * pow(dot(-vRay, refl), object.shinniness);

        if (object.kind == OBJ_KIND_PLANE) {
            diffuse = abs(diffuse);
        }
        diffuse = clampColor(diffuse);
        specular = clampColor(specular);
        color += specular;
        color += diffuse;
    }

    return clampColor(color);
}

#define REFRACTION_INDEX_NORMAL 1
#define REFRACTION_INDEX_SPHERE 1.5
#define MAX_TRACING_COUNT 5
void bounceLightRay(inout StraightLine ray, out Intersection intersection) {
    int i;
    float refractionIndex = REFRACTION_INDEX_NORMAL;
    intersection = findFirstIntersectingObject(ray.p0, ray.v, true, -1, false);
    for (i = 0; i < MAX_TRACING_COUNT; i++) {
        if (intersection.objectIndex < 0) {
            break;
        }
        
        if ((getObjectFlags(intersection.objectIndex) & OBJ_FLAGS_REFLECTIVE) != 0) {
            vec3 vRay = normalize(reflect(ray.v, intersection.pointNormal));
            ray = StraightLine(intersection.point, vRay);
        }
        else if ((getObjectFlags(intersection.objectIndex) & OBJ_FLAGS_TRANSPARENT) != 0) {
            vec4 object = objects[intersection.objectIndex];
            // Refract only in case of a transparent sphere, but not if we are tanget to the radius
            if (isObjectOfKind(object, OBJ_KIND_SPHERE) && !isZero(dot(ray.v, intersection.pointNormal))) {
                // TODO: calculate the next refraction index based on the object hit:
                //   * if we hit the same shpere we were before (but now from the inside), set to 1.
                //   * otherwise, based on whether the object we hit is transparent or not.
                // TODO: calculate the initial refraction (the very first) index based on the object we start in (if we do)
                float nextRefractionIndex = (REFRACTION_INDEX_NORMAL + REFRACTION_INDEX_SPHERE) - refractionIndex;
                float refractionRatio = refractionIndex / nextRefractionIndex;
                vec3 vRay = normalize(refract(ray.v, intersection.pointNormal, refractionRatio));
                ray = StraightLine(intersection.point, vRay);
                refractionIndex = nextRefractionIndex;
            }
            else {
                // This is a thin surface, so the light goes through it without changes.
                // Therefore, this shouldn't count, so accomodate here for the increment that
                // comes after this.
                i--;
            }
        }
        else {
            break;
        }

        intersection = findFirstIntersectingObject(ray.p0, ray.v, true, intersection.objectIndex, false);
    }
}

void main()
{
    vec3 vRay = normalize(position0.xyz - eye.xyz);
    StraightLine ray = StraightLine(eye.xyz, vRay);
    Intersection intersection;
    bounceLightRay(ray, intersection);

    vec3 color = vec3(0);
    if (intersection.objectIndex >= 0) {
        color = calculateColor_noTracing(ray.v, intersection);
    }
    outColor = vec4(color, 1);
}
