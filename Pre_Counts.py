import os
import shutil

# specify the directory you want to use
directory = r'C:\Users\ADMIN\Desktop\Unripe_ Assembled transcripts'

for filename in os.listdir(directory):
    if os.path.isfile(os.path.join(directory, filename)):
        # get the file name without the extension
        folder_name = os.path.splitext(filename)[0]
        # create a new folder for each file
        new_folder_path = os.path.join(directory, folder_name)
        os.makedirs(new_folder_path, exist_ok=True)
        
        # move the file to the new folder
        shutil.move(os.path.join(directory, filename), os.path.join(new_folder_path, filename))
