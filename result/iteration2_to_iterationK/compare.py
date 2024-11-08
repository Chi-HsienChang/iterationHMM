import matplotlib.pyplot as plt
import matplotlib.image as mpimg

# List of image indices
indices = range(1, 7)

# Create a figure with 6 rows and 2 columns
fig, axes = plt.subplots(nrows=6, ncols=2, figsize=(10, 30))
# Add titles to the top subplots
axes[0, 0].set_title('wo removing insertion', fontsize=16)
axes[0, 1].set_title('w removing insertion', fontsize=16)

# Adjust the spacing between subplots
plt.subplots_adjust(hspace=0.01, wspace=0.01)

for idx, i in enumerate(indices):
    # Load images
    img_w = mpimg.imread(f'w_CCCHCCHH_i{i}.png')
    img_wo = mpimg.imread(f'i{i}.png')
    
    # Display w_iX.png in the left column
    axes[idx, 0].imshow(img_w)
    axes[idx, 0].axis('off')  # Hide axis ticks
    
    # Display wo_iX.png in the right column
    axes[idx, 1].imshow(img_wo)
    axes[idx, 1].axis('off')  # Hide axis ticks

# Optionally, save the combined figure
plt.savefig('combined_figure_1108.png', bbox_inches='tight', dpi=300)
# plt.show()
