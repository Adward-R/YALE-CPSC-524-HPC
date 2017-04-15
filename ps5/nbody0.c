#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

struct BodyType { 
  float x, y, z;
  float vx, vy, vz; 
};

void MoveBodies(const int nBodies, struct BodyType* const bodies, const float dt) {

  // Avoid singularity and interaction with self
  const float softening = 1e-20;

  // Loop over bodies that experience force
  for (int i = 0; i < nBodies; i++) { 

    // Components of the gravity force on body i
    float Fx = 0, Fy = 0, Fz = 0; 
      
    // Loop over bodies that exert force: vectorization expected here
    for (int j = 0; j < nBodies; j++) { 
      
      // Newton's law of universal gravity
      const float dx = bodies[j].x - bodies[i].x;
      const float dy = bodies[j].y - bodies[i].y;
      const float dz = bodies[j].z - bodies[i].z;
      const float drSquared  = dx*dx + dy*dy + dz*dz + softening;
      const float drPower32  = pow(drSquared, 3.0/2.0);
	
      // Calculate the net force
      Fx += dx / drPower32;  
      Fy += dy / drPower32;  
      Fz += dz / drPower32;

    }

    // Accelerate bodies in response to the gravitational force
    bodies[i].vx += dt*Fx; 
    bodies[i].vy += dt*Fy; 
    bodies[i].vz += dt*Fz;
  }

  // Move bodies according to their velocities
  for (int i = 0 ; i < nBodies; i++) { 
    bodies[i].x  += bodies[i].vx*dt;
    bodies[i].y  += bodies[i].vy*dt;
    bodies[i].z  += bodies[i].vz*dt;
  }
}

int main(const int argc, const char** argv) {

  // Problem size and other parameters
  const int nBodies = (argc > 1 ? atoi(argv[1]) : 16384);
  const int nSteps = 10;  // Duration of test
  const float dt = 0.01f; // Body propagation time step

  // Body data stored as an Array of Structures (AoS)
  struct BodyType bodies[nBodies];

  // Initialize random number generator and bodies
  srand(0);
  float randmax;
  randmax = (float) RAND_MAX;
  for(int i = 0; i < nBodies; i++) {
    bodies[i].x = ((float) rand())/randmax; 
    bodies[i].y = ((float) rand())/randmax; 
    bodies[i].z = ((float) rand())/randmax; 
    bodies[i].vx = ((float) rand())/randmax; 
    bodies[i].vy = ((float) rand())/randmax; 
    bodies[i].vz = ((float) rand())/randmax; 
  }

  // Compute initial center of mass  
  float comx = 0.0f, comy=0.0f, comz=0.0f;
  for (int i=0; i<nBodies; i++) {
    comx += bodies[i].x;
    comy += bodies[i].y;
    comz += bodies[i].z;
  }
  comx = comx / nBodies;
  comy = comy / nBodies;
  comz = comz / nBodies;

  printf("Initial center of mass: (%g, %g, %g)\n", comx, comy, comz);

  // Perform benchmark
  printf("\n\033[1mNBODY Version 00\033[0m\n");
  printf("\nPropagating %d bodies using 1 thread on %s...\n\n", 
	 nBodies, "CPU");

  double rate = 0, dRate = 0; // Benchmarking data
  const int skipSteps = 3; // Set this to a positive int to skip warm-up steps

  printf("\033[1m%5s %10s %10s %8s\033[0m\n", "Step", "Time, s", "Interact/s", "GFLOP/s"); fflush(stdout);

  for (int step = 1; step <= nSteps; step++) {

    const double tStart = omp_get_wtime(); // Start timing
    MoveBodies(nBodies, bodies, dt);
    const double tEnd = omp_get_wtime(); // End timing

    // These are for calculating flop rate. It ignores symmetry and 
    // estimates 20 flops per body-body interaction in MoveBodies
    const float HztoInts   = (float)nBodies * (float)(nBodies-1) ;
    const float HztoGFLOPs = 20.0*1e-9*(float)nBodies*(float)(nBodies-1);

    if (step > skipSteps) { 
      // Collect statistics 
      rate  += HztoGFLOPs / (tEnd - tStart); 
      dRate += HztoGFLOPs * HztoGFLOPs / ((tEnd-tStart)*(tEnd-tStart)); 
    }

    printf("%5d %10.3e %10.3e %8.1f %s\n", 
	   step, (tEnd-tStart), HztoInts/(tEnd-tStart), HztoGFLOPs/(tEnd-tStart), (step<=skipSteps?"*":""));
    fflush(stdout);
  }

  rate/=(double)(nSteps-skipSteps); 
  dRate=sqrt(fabs(dRate/(double)(nSteps-skipSteps)-rate*rate));

  printf("-----------------------------------------------------\n");
  printf("\033[1m%s %4s \033[42m%10.1f +- %.1f GFLOP/s\033[0m\n",
	 "Average performance:", "", rate, dRate);
  printf("-----------------------------------------------------\n");
  printf("* - warm-up, not included in average\n\n");

  // Compute final center of mass
  comx = 0.0f; comy=0.0f; comz=0.0f;
  for (int i=0; i<nBodies; i++) {
    comx += bodies[i].x;
    comy += bodies[i].y;
    comz += bodies[i].z;
  }
  comx = comx / nBodies;
  comy = comy / nBodies;
  comz = comz / nBodies;
  
  printf("Final center of mass: (%g, %g, %g)\n", comx, comy, comz);

}
