from PIL import Image
import glob
import math
import os

# Use the current working directory for the folder path
folder_path = '.'

# Get a list of all PNG image files in the current directory
image_files = glob.glob(f'{folder_path}/*.png')

# Dynamically calculate the number of rows and columns for the grid
n = len(image_files)
cols = int(math.ceil(math.sqrt(n)))  # Square root to get a balanced grid
rows = int(math.ceil(n / cols))

# Assuming all images are the same size, open the first image to get its size
sample_image = Image.open(image_files[0])
width, height = sample_image.size

# Create a new blank image canvas large enough to contain all subplots
canvas_width = cols * width
canvas_height = rows * height
combined_image = Image.new('RGB', (canvas_width, canvas_height))

# Paste each image into its correct position on the canvas
for index, image_file in enumerate(image_files):
    img = Image.open(image_file)
    # Calculate the position based on the index
    x_offset = (index % cols) * width
    y_offset = (index // cols) * height
    combined_image.paste(img, (x_offset, y_offset))

# Save the combined image
output_path = os.path.join(folder_path, 'combined_venn_diagrams.png')
combined_image.save(output_path)
print(f"Combined image saved to {output_path}")

# Optionally, display the combined image
# combined_image.show()
