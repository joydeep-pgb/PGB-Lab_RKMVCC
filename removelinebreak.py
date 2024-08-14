# Open the original text file for reading
with open('GF_RF.interaction.txt', 'r') as original_file:
    # Read the entire content of the original text file
    original_content = original_file.read()

# Remove line breaks from the content
content_without_line_breaks = original_content.replace('\n', '')

# Open a new file for writing the content without line breaks
with open('GF_RF.interaction.txt', 'w') as arranged_file:
    # Write the content without line breaks to the new file
    arranged_file.write(content_without_line_breaks)
