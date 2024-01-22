# This script provides some of the numbers written in the manuscript

# INPUT: none
# OUTPUT: none

# >>>>>>>>>

# Calculate effective slope if the absorption slope is 0.54 and 1/3 of the total ‘vital effect’ is induced by diffusion
effective_slope = 2/3 * 0.54 + 0.505 * 1/3
print(f"The vital effect slope is {effective_slope:.3f} if 1/3 of the total ‘vital effect’ is induced by diffusion")

# Calculate diffusion-induced 'vital effect' percentage
coral_slope = 0.530
abs_values = [0.538, 0.541]
for abs in abs_values:
    x = (1 - ((coral_slope - 0.505) / (abs - 0.505))) * 100
    print(f"\nIf the absorption slope is {abs}, then {x:.0f}% of the total ‘vital effect’ is induced by diffusion")


# Calculate absorption slopes with revised theta estimates from Bajnai et al. (2023)
hydrox_contrib = [0.5, 0.8]
abs_values = []
for hydrx in hydrox_contrib:
    abs = 0.532 * hydrx + 0.526 * (1-hydrx)
    abs_values.append(abs)  # Append the calculated abs value to the list
    print(f"\nThe absorption slope is {abs:.3f}, with {hydrx*100:.0f}% of the total ‘vital effect’ induced by hydroxylation")

# Calculate diffusion-induced 'vital effect' percentage with revised absorption theta estimates from Bajnai et al. (2023)
coral_slopes = [0.528, 0.532]
diff_percentages = []
for abs in abs_values:
    for coral_slope in coral_slopes:
        x = (1 - ((coral_slope - 0.505) / (abs - 0.505))) * 100
        diff_percentages.append(x)

print(f"\nUp to {max(diff_percentages):.0f}% of the total ‘vital effect’ can be induced by diffusion")