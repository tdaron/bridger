#version 400 core
layout (location = 0) in vec3 aPos;
// layout (location = 1) in vec3 color; // Unused
layout (location = 2) in vec2 displacement;

out vec3 vertexColor; // Output color to fragment shader

// Scaling factor for the displacement effect on position
const float FACTOR = 3000.0f;

// --- TUNABLE PARAMETER ---
// Represents the scaled displacement magnitude that should result in full red.
// Adjust as needed for your data range.
const float MAX_EXPECTED_SCALED_MAGNITUDE = 0.01f;
// -------------------------

// Define key colors for the new gradient
const vec3 COLOR_ZERO = vec3(0.0, 0.0, 1.0); // Blue
const vec3 COLOR_LOW  = vec3(0.0, 1.0, 1.0); // Cyan
const vec3 COLOR_MID  = vec3(1.0, 1.0, 0.0); // Yellow
const vec3 COLOR_HIGH = vec3(1.0, 0.0, 0.0); // Red

// Define thresholds for color transitions (approx 0.0 -> 1/3 -> 2/3 -> 1.0)
const float THRESHOLD1 = 1.0 / 3.0;
const float THRESHOLD2 = 2.0 / 3.0;

void main()
{
    // Calculate the displaced vertex position
    vec3 displacedPos = vec3(aPos.x + displacement.x * FACTOR, aPos.y + displacement.y * FACTOR, aPos.z);
    gl_Position = vec4(displacedPos, 1.0);

    // --- Color Calculation ---

    // Calculate the magnitude (length) of the scaled displacement vector
    float displacementMagnitude = length(displacement * FACTOR);

    // Normalize the magnitude to a [0.0, 1.0] range.
    float normMag = clamp(displacementMagnitude / MAX_EXPECTED_SCALED_MAGNITUDE, 0.0, 1.0);

    // Interpolate Blue -> Cyan -> Yellow -> Red
    if (normMag < THRESHOLD1)
    {
        // Segment 1: Blue to Cyan
        // Remap normMag from [0, THRESHOLD1) to [0, 1)
        // factor = normMag / THRESHOLD1 which is normMag * 3.0
        float factor = normMag * 3.0;
        vertexColor = mix(COLOR_ZERO, COLOR_LOW, factor);
    }
    else if (normMag < THRESHOLD2)
    {
        // Segment 2: Cyan to Yellow
        // Remap normMag from [THRESHOLD1, THRESHOLD2) to [0, 1)
        // factor = (normMag - THRESHOLD1) / (THRESHOLD2 - THRESHOLD1)
        // factor = (normMag - 1/3) / (1/3) which is (normMag - 1/3) * 3.0
        float factor = (normMag - THRESHOLD1) * 3.0;
        vertexColor = mix(COLOR_LOW, COLOR_MID, factor);
    }
    else
    {
        // Segment 3: Yellow to Red
        // Remap normMag from [THRESHOLD2, 1.0] to [0, 1]
        // factor = (normMag - THRESHOLD2) / (1.0 - THRESHOLD2)
        // factor = (normMag - 2/3) / (1/3) which is (normMag - 2/3) * 3.0
        float factor = (normMag - THRESHOLD2) * 3.0;
        vertexColor = mix(COLOR_MID, COLOR_HIGH, factor);
    }
}
