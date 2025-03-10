#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <log.h>
#include "mesh.h"

Mesh* readMesh(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        log_error("Error opening file: %s", filename);
        exit(1);
        return NULL;
    }
    
    Mesh* mesh = (Mesh*)malloc(sizeof(Mesh));
    if (!mesh) {
        fclose(file);
        return NULL;
    }
    
    Vertex* vertices = NULL;
    Triangle* triangles = NULL;
    triangles = NULL;
    mesh->vertexArray = NULL;
    
    char line[256];
    
    // Read number of nodes
    if (fgets(line, sizeof(line), file) == NULL) {
        log_error("Error reading number of nodes");
        free(mesh);
        fclose(file);
        return NULL;
    }
    
    int numNodes;
    if (sscanf(line, "Number of nodes %d", &numNodes) != 1) {
        log_error("Invalid format for number of nodes");
        free(mesh);
        fclose(file);
        return NULL;
    }
    
    log_info("Reading %d nodes...", numNodes);
    
    mesh->numVertices = numNodes;
    vertices = (Vertex*)malloc(numNodes * sizeof(Vertex));
    if (!vertices) {
        log_error("Memory allocation error for vertices");
        free(mesh);
        fclose(file);
        return NULL;
    }
    
    // Read node coordinates
    for (int i = 0; i < numNodes; i++) {
        if (fgets(line, sizeof(line), file) == NULL) {
            log_error("Error reading node %d", i);
            free(vertices);
            free(mesh);
            fclose(file);
            return NULL;
        }
        
        int nodeIdx;
        float x, y;
        if (sscanf(line, "%d : %e %e", &nodeIdx, &x, &y) != 3) {
            log_error("Invalid format for node %d", i);
            free(vertices);
            free(mesh);
            fclose(file);
            return NULL;
        }
        
        // Store node data
        vertices[i].x = x;
        vertices[i].y = y;
        vertices[i].z = 0.0f;  // Set z to 0 for 2D mesh
    }
    
    // Skip to triangles section
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Number of triangles")) {
            break;
        }
    }
    
    // Read number of triangles
    int numTriangles;
    if (sscanf(line, "Number of triangles %d", &numTriangles) != 1) {
        log_error("Invalid format for number of triangles");
        free(vertices);
        free(mesh);
        fclose(file);
        return NULL;
    }
    
    log_info("Reading %d triangles...", numTriangles);
    
    mesh->numTriangles = numTriangles;
    triangles = (Triangle*)malloc(numTriangles * sizeof(Triangle));
    if (!triangles) {
        log_error("Memory allocation error for triangles");
        free(vertices);
        free(mesh);
        fclose(file);
        return NULL;
    }
    
    // Read triangle vertex indices
    int triCount = 0;
    while (fgets(line, sizeof(line), file) && triCount < numTriangles) {
        int triIdx, v1, v2, v3;
        if (sscanf(line, "%d : %d %d %d", &triIdx, &v1, &v2, &v3) != 4) {
            // Skip if not a triangle line (might be header or other content)
            continue;
        }
        
        // Store triangle data with 0-based indexing
        triangles[triCount].v1 = v1;
        triangles[triCount].v2 = v2;
        triangles[triCount].v3 = v3;
        triCount++;
    }
    
    if (triCount != numTriangles) {
        log_error("Warning: Read only %d triangles out of %d", triCount, numTriangles);
        mesh->numTriangles = triCount;
    }
    
    fclose(file);
    
    log_info("Normalizing mesh coordinates...");
    
    // Find mesh bounding box
    float minX = FLT_MAX, maxX = -FLT_MAX;
    float minY = FLT_MAX, maxY = -FLT_MAX;
    
    for (int i = 0; i < mesh->numVertices; i++) {
        if (vertices[i].x < minX) minX = vertices[i].x;
        if (vertices[i].x > maxX) maxX = vertices[i].x;
        if (vertices[i].y < minY) minY = vertices[i].y;
        if (vertices[i].y > maxY) maxY = vertices[i].y;
    }
    
    // Calculate normalization parameters
    float width = maxX - minX;
    float height = maxY - minY;
    float scale = 1.9f / (width > height ? width : height); // Scale to fit in [-0.95, 0.95]
    float centerX = (minX + maxX) / 2.0f;
    float centerY = (minY + maxY) / 2.0f;
    
    log_info("Mesh bounds: X[%.2f, %.2f], Y[%.2f, %.2f]", minX, maxX, minY, maxY);
    log_info("Center: (%.2f, %.2f), Scale: %.5f", centerX, centerY, scale);
    
    // Create vertex array for OpenGL (3 vertices per triangle, 3 coordinates per vertex)
    mesh->vertexArraySize = mesh->numTriangles * 3 * 3;
    mesh->vertexArray = (float*)malloc(mesh->vertexArraySize * sizeof(float));
    if (!mesh->vertexArray) {
        log_error("Memory allocation error for vertex array");
        free(triangles);
        free(vertices);
        free(mesh);
        return NULL;
    }
    
    // Fill the vertex array with normalized coordinates
    for (int i = 0; i < mesh->numTriangles; i++) {
        Triangle* tri = &triangles[i];
        
        // First vertex
        mesh->vertexArray[i*9 + 0] = (vertices[tri->v1].x - centerX) * scale;
        mesh->vertexArray[i*9 + 1] = (vertices[tri->v1].y - centerY) * scale;
        mesh->vertexArray[i*9 + 2] = vertices[tri->v1].z;
        
        // Second vertex
        mesh->vertexArray[i*9 + 3] = (vertices[tri->v2].x - centerX) * scale;
        mesh->vertexArray[i*9 + 4] = (vertices[tri->v2].y - centerY) * scale;
        mesh->vertexArray[i*9 + 5] = vertices[tri->v2].z;
        
        // Third vertex
        mesh->vertexArray[i*9 + 6] = (vertices[tri->v3].x - centerX) * scale;
        mesh->vertexArray[i*9 + 7] = (vertices[tri->v3].y - centerY) * scale;
        mesh->vertexArray[i*9 + 8] = vertices[tri->v3].z;
    }
    
    // Print debug info for first triangle
    log_info("First triangle vertices after normalization:");
    log_info("V1: (%.2f, %.2f, %.2f)", 
           mesh->vertexArray[0], mesh->vertexArray[1], mesh->vertexArray[2]);
    log_info("V2: (%.2f, %.2f, %.2f)", 
           mesh->vertexArray[3], mesh->vertexArray[4], mesh->vertexArray[5]);
    log_info("V3: (%.2f, %.2f, %.2f)", 
           mesh->vertexArray[6], mesh->vertexArray[7], mesh->vertexArray[8]);
    
    log_info("Mesh loaded successfully: %d vertices, %d triangles", 
           mesh->numVertices, mesh->numTriangles);
    
    free(triangles);
    free(vertices);
    return mesh;
}

// Function to free mesh resources
void freeMesh(Mesh* mesh) {
    if (mesh) {
        if (mesh->vertexArray) free(mesh->vertexArray);
        free(mesh);
    }
}
