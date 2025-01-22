from PIL import Image
import os

# List of image filenames
image_files = [f"i{n}.png" for n in range(1, 11)]

# Open images and store them in a list
images = [Image.open(img) for img in image_files]

# Determine size of a single image
width, height = images[0].size

# Option to stack images horizontally or vertically
mode = 'horizontal'  # Change to 'vertical' if needed

if mode == 'horizontal':
    # Create a new image with a width equal to width of all images combined
    total_width = width * len(images)
    new_image = Image.new('RGB', (total_width, height))
    x_offset = 0
    for img in images:
        new_image.paste(img, (x_offset, 0))
        x_offset += img.width
else:
    # Create a new image with a height equal to height of all images combined
    total_height = height * len(images)
    new_image = Image.new('RGB', (width, total_height))
    y_offset = 0
    for img in images:
        new_image.paste(img, (0, y_offset))
        y_offset += img.height

# Save the merged image
new_image.save('merged_image.png')
print("Images have been merged and saved as 'merged_image.png'")
