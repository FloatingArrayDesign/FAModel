
import numpy as np
import matplotlib.pyplot as plt
from installation_dynamic import getInstallationDynamic
from scipy.stats import lognorm

def getInstallationDynamicMC(profile_map, location, D1, D2, L1, L2, ballast, drop_height, N_sim=2000):
    np.random.seed(16)
    
    # Lognormal distribution parameters
    Suk_mean, Suk_std = 1.9, 0.9
    Sti_mean, Sti_std = 3.2, 1.0

    Suk_samples = np.random.lognormal(mean=np.log(Suk_mean), sigma=Suk_std/Suk_mean, size=N_sim)
    Sti_samples = np.random.lognormal(mean=np.log(Sti_mean), sigma=Sti_std/Sti_mean, size=N_sim)

    final_depths = []

    for Suk, Sti in zip(Suk_samples, Sti_samples):
        # Update profile with consistent linear Su(z)
        profile = [dict(profile_map[0])]  
        Su0 = profile[0]['layers'][0]['Su_top']
        for layer in profile[0]['layers']:
            z_top = layer['top']
            z_bot = layer['bottom']
            layer['Su_top'] = Su0 + Suk*z_top
            layer['Su_bot'] = Su0 + Suk*z_bot
        
        # Override getInstallationDynamic with fixed Sti
        result = getInstallationDynamic(profile, location, D1, D2, L1, L2, ballast, drop_height, plot=False)
        
        final_depths.append(result['final_depth'])
        
    # Fit lognormal distribution to the data
    shape, loc, scale = lognorm.fit(final_depths, floc=0)  
    x = np.linspace(min(final_depths), max(final_depths), 500)
    pdf = lognorm.pdf(x, shape, loc=loc, scale=scale)

    # Create subplot
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    
    # Fit lognormals
    s_Suk, loc_Suk, scale_Suk = lognorm.fit(Suk_samples, floc=0)
    s_Sti, loc_Sti, scale_Sti = lognorm.fit(Sti_samples, floc=0)
        
    # Generate x values for plotting
    x_Suk = np.linspace(min(Suk_samples), max(Suk_samples), 500)
    x_Sti = np.linspace(min(Sti_samples), max(Sti_samples), 500)
    
    # Compute PDFs
    pdf_Suk = lognorm.pdf(x_Suk, s_Suk, loc=loc_Suk, scale=scale_Suk)
    pdf_Sti = lognorm.pdf(x_Sti, s_Sti, loc=loc_Sti, scale=scale_Sti)
    
    # Compute max density from histograms and PDFs for setting common ylim
    hist_Suk_vals, _ = np.histogram(Suk_samples, bins=40, density=True)
    hist_Sti_vals, _ = np.histogram(Sti_samples, bins=40, density=True)
    y_max = max(max(hist_Suk_vals), max(hist_Sti_vals), max(pdf_Suk), max(pdf_Sti)) * 1.1
    
    # Plot Suk PDF
    axs[0].hist(Suk_samples, bins=40, density=True, alpha=0.6, color='lightgreen', edgecolor='g', label='Samples')
    axs[0].plot(x_Suk, pdf_Suk, 'r', label='Fitted Lognormal')
    axs[0].set_title('Lognormal Fit for $S_{uk}$')
    axs[0].set_xlabel('Suk (kPa/m)')
    axs[0].set_ylabel('Density')
    axs[0].set_ylim(0, y_max)
    axs[0].legend()
    axs[0].grid(True)
    
    # Plot Sti PDF
    axs[1].hist(Sti_samples, bins=40, density=True, alpha=0.6, color='lightblue', edgecolor='b', label='Samples')
    axs[1].plot(x_Sti, pdf_Sti, 'r', label='Fitted Lognormal')
    axs[1].set_title('Lognormal Fit for $S_{ti}$')
    axs[1].set_xlabel('Sti (dimensionless)')
    axs[1].set_ylabel('Density')
    axs[1].set_ylim(0, y_max)
    axs[1].legend()
    axs[1].grid(True)
       
    # Plot histogram with fitted lognormal
    plt.figure(figsize=(8, 5))
    plt.hist(final_depths, bins=40, alpha=0.6, density=True, color='lightsalmon', edgecolor='r', label='Simulated PDF')
    plt.plot(x, pdf, 'r', lw=2, label='Fitted Lognormal')
    plt.axvline(np.median(final_depths), color='k', linestyle='--', label=f'Median = {np.median(final_depths):.2f} m')
    plt.title('Monte Carlo Simulation with Fitted Lognormal Distribution')
    plt.xlabel('Penetration (m)')
    plt.ylabel('Density')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Print parameters
    print(f'Lognormal fit parameters:')
    print(f'  mean                   = {np.mean(final_depths):.2f}')
    print(f'  mean - 1 std deviation = {np.mean(final_depths) - np.std(final_depths):.2f}')
    print(f'  mean + 1 std deviation = {np.mean(final_depths) + np.std(final_depths):.2f}')

    return final_depths

# Example usage in main
if __name__ == '__main__':
    profile_map = [
        {
            'name': 'CPT_1',
            'x': 0, 'y': 0,
            'layers': [
                {'top': 0.0, 'bottom': 30.0, 'soil_type': 'clay',
                 'gamma_top': 8.0, 'gamma_bot': 8.5,
                 'Su_top': 5, 'Su_bot': 20},
                {'top': 30.0, 'bottom': 180.0, 'soil_type': 'clay',
                 'gamma_top': 8.5, 'gamma_bot': 9.0,
                 'Su_top': 20, 'Su_bot': 150}
            ]
        }
    ]

    location = 'CPT_1'
    D1 = 3.0
    D2 = 1.2
    L1 = 5.0
    L2 = 5.0
    ballast = 300000
    drop_height = 200

    # Run Monte Carlo Simulation
    depths = getInstallationDynamicMC(profile_map, location, D1, D2, L1, L2, ballast, drop_height)
