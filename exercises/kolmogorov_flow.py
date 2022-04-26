# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 16:25:19 2022

Reimplemented by Joshua Brinkerhoff (@okcfdlab) into Python from a Julia code developed by Felix Köhler (@Ceyron):
https://github.com/Ceyron/machine-learning-and-simulation/blob/main/english/simulation_scripts/kolmogorov_turbulence.jl

Solves the equations of fluid flow using a "Stable Fluids"-like algorithm
based on the FFT in a Kolmogorov Flow scenario. This setup results in forced
isotropic turbulence.
Momentum:           ∂u/∂t + (u ⋅ ∇) u = − ∇p + ν ∇²u + f
Incompressibility:  ∇ ⋅ u = 0
u:  Velocity (2d vector)
p:  Pressure
f:  Forcing (here Kolmogorov Forcing)
ν:  Kinematic Viscosity
t:  Time
∇:  Nabla operator (defining nonlinear convection, gradient and divergence)
∇²: Laplace Operator
----
Setup:
A rectangular domain (16x9) with periodic boundary conditions in both axes.

-> The fluid is initialized at rest
-> Then, a forcing is applied according to
        ⎡   c * sin(k * y)  ⎤
    f = ⎢                   ⎥
        ⎣         0         ⎦
        with scaling c and wavenumber k
----- 
Solution Strategy:
-> Start with zero velocity everywhere: u = [0, 0]
1. Add forces
    w₁ = u + Δt f
2. Convect by self-advection (set the value at the current
   location to be the value at the position backtraced
   on the streamline.) -> unconditionally stable
    w₂ = w₁(p(x, −Δt))
3. Subtract mean velocity for stabilization
4. Diffuse and Project in Fourier Domain
    4.1 Forward Transformation into Fourier Domain
        w₂ → w₃
    
    4.2 Diffuse by "low-pass filtering" (convolution
        is multiplication in the Fourier Domain)
        w₄ = exp(− k² ν Δt) w₃
    
    4.3 Compute the (pseudo-) pressure in the Fourier Domain
        by evaluating the divergence in the Fourier Domain
        q = w₄ ⋅ k / ||k||₂
    
    4.4 Correct the velocities such that they are incompressible
        w₅ = w₄ − q k / ||k||₂
    
    4.5 Inverse Transformation back into spatial domain
        w₆ ← w₅
5. Subtract mean velocity for stabilization
6. Repeat
k = [ k_x, k_y ] are the spatial frequencies (= wavenumbers)
The Fourier Transformation implicitly prescribes the periodic
Boundary Conditions
-----
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator


N_POINTS_Y = 360
ASPECT_RATIO = 16.0/9.0
KINEMATIC_VISCOSITY = 1.0 / 1000.0
TIME_STEP_LENGTH = 0.001
N_TIME_STEPS = 750

FORCING_WAVENUMBER = 8
FORCING_SCALE = 100.0


n_points_x = int(N_POINTS_Y*ASPECT_RATIO)
x_extent = n_points_x/N_POINTS_Y - 1.0e-7
x_interval = np.linspace(0.0,x_extent,n_points_x)
y_interval = np.linspace(0.0,1.0,N_POINTS_Y)

coords_x, coords_y = np.meshgrid(x_interval, y_interval, indexing='ij')

wavenumbers_1d_x = np.fft.fftfreq(n_points_x) * n_points_x
n_fft_points_x = np.size(wavenumbers_1d_x)
wavenumbers_1d_y = np.fft.rfftfreq(N_POINTS_Y) * N_POINTS_Y
n_fft_points_y = np.size(wavenumbers_1d_y)

wavenumbers_x, wavenumbers_y = np.meshgrid(wavenumbers_1d_x, wavenumbers_1d_y, indexing='ij')

# compute matrix of wavenumber norms and convert zero elements in wavenumber_norm matrix to 1.0
wavenumbers_norm = np.copy(wavenumbers_x)
for kx in range(n_fft_points_x):
    for ky in range(n_fft_points_y):
        wavenumbers_norm[kx,ky] = np.linalg.norm([wavenumbers_1d_x[kx],wavenumbers_1d_y[ky]])
        if wavenumbers_norm[kx,ky] == 0:
            wavenumbers_norm[kx,ky] = 1

decay = np.exp(-TIME_STEP_LENGTH * KINEMATIC_VISCOSITY * wavenumbers_norm**2)


normalized_wavenumbers_x = wavenumbers_x/wavenumbers_norm
normalized_wavenumbers_y = wavenumbers_y/wavenumbers_norm

force_x = FORCING_SCALE * np.sin(FORCING_WAVENUMBER * np.pi * coords_y)

# Preallocate arrays
backtraced_coordinates_x = np.zeros([n_points_x, N_POINTS_Y])
backtraced_coordinates_y = np.zeros([n_points_x, N_POINTS_Y])

velocity_x = np.zeros([n_points_x, N_POINTS_Y])
velocity_y = np.zeros([n_points_x, N_POINTS_Y])

velocity_x_prev = np.zeros([n_points_x, N_POINTS_Y])
velocity_y_prev = np.zeros([n_points_x, N_POINTS_Y])

velocity_x_fft = np.zeros([n_fft_points_x, n_fft_points_y], dtype=complex)
velocity_y_fft = np.zeros([n_fft_points_x, n_fft_points_y], dtype=complex)
pressure_fft = np.zeros([n_fft_points_x, n_fft_points_y], dtype=complex)

dudy_fft =  np.zeros([n_fft_points_x, n_fft_points_y], dtype=complex)
dvdx_fft =  np.zeros([n_fft_points_x, n_fft_points_y], dtype=complex)
curl_fft =  np.zeros([n_fft_points_x, n_fft_points_y], dtype=complex)
curl = np.zeros([n_points_x, N_POINTS_Y])


# Timestepping loop

for iter in range(N_TIME_STEPS):
    
    if (iter % 5) == 0:
        print("Timestep: " + str(iter))
        
    # (1) Apply the forces
    velocity_x_prev += force_x * TIME_STEP_LENGTH

    # (2) Self-advection by backtracing and interpolation
    backtraced_coordinates_x = np.mod(
        coords_x - TIME_STEP_LENGTH*velocity_x_prev,
        x_extent
    )

    backtraced_coordinates_y = np.mod(
        coords_y - TIME_STEP_LENGTH*velocity_y_prev,
        1.0
    )

    x_interpolator = RegularGridInterpolator((x_interval, y_interval), velocity_x_prev)
    y_interpolator = RegularGridInterpolator((x_interval, y_interval), velocity_y_prev)
    backtraced_locs = np.column_stack((backtraced_coordinates_x.flatten(),backtraced_coordinates_y.flatten()))
   
    # Interpolate velocities at backtraced locations. Note, the output is a flattened (1D) array
    velocity_x = x_interpolator(backtraced_locs)
    velocity_y = y_interpolator(backtraced_locs)
    

    # (3) Stabilize by subtracting mean velocity    
    velocity_x -= np.mean(velocity_x)
    velocity_y -= np.mean(velocity_y)
    
    # Resize velocity arrays
    velocity_x.resize(n_points_x,N_POINTS_Y)
    velocity_y.resize(n_points_x,N_POINTS_Y)

    # (4.1) Transform into wavenumber space
    velocity_x_fft = np.fft.rfft2(velocity_x)
    velocity_y_fft = np.fft.rfft2(velocity_y)
    
    # (4.2) Diffuse by low-pass filtering
    velocity_x_fft *= decay
    velocity_y_fft *= decay
    
    # (4.3) Compute the pseudo pressure via divergence. In the wavenumber space, divergence transforms
    #       to multiplication of the normalized wavenumbers
    pressure_fft = (
        velocity_x_fft * normalized_wavenumbers_x
        +
        velocity_y_fft * normalized_wavenumbers_y
        )
    
    # (4.4) Project the velocity to be incompressible
    velocity_x_fft -= pressure_fft * normalized_wavenumbers_x
    velocity_y_fft -= pressure_fft * normalized_wavenumbers_y
    
    # (4.5) Transform back into the spatial domain
    velocity_x = np.fft.irfft2(velocity_x_fft, velocity_x.shape)
    velocity_y = np.fft.irfft2(velocity_y_fft, velocity_y.shape)
    
    # (5) Stabilize by subtracting the mean velocities
    velocity_x -= np.mean(velocity_x)
    velocity_y -= np.mean(velocity_y)
    
    # (6) Advance in time
    velocity_x_prev = np.copy(velocity_x)
    velocity_y_prev = np.copy(velocity_y)
    
    # (7) Compute the velocity curl
    dudy_fft = wavenumbers_y * velocity_x_fft * 1j
    dvdx_fft = wavenumbers_x * velocity_y_fft * 1j
    curl_fft = dvdx_fft - dudy_fft
    curl = np.fft.irfft2(curl_fft, curl.shape)
    
    curl = np.sign(curl) * np.sqrt(np.abs(curl) / np.quantile(curl.flatten(), 0.8))
    
    fig, ax = plt.subplots()
    cnt = ax.contourf(
        coords_x,
        coords_y,
        curl,
        cmap='YlGnBu',
        levels=100
        )
    ax.axis('off')
    labeltext = 'Iter: ' + str(iter).zfill(5)
    ax.text(0.1, 0.9, labeltext, fontsize=12)
    fig.set_size_inches(5*ASPECT_RATIO,5)

    iter_string = str(iter)
    filename = 'kolmogorov_' + iter_string.zfill(5) + '.png'
    fig.savefig(filename, bbox_inches='tight', transparent=True, dpi=96)
    plt.close()
    plt.clf()
    
