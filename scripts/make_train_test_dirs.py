import os
import shutil
import random
import numpy as np

path_name = '../test_gaussian_parallel/'
file_names = os.listdir(path_name)


# shuffle and split data
np.random.shuffle(file_names)
train_dirs = file_names[:int(len(file_names) * 0.8)]
test_dirs = file_names[int(len(file_names) * 0.8):]


train_copy_path = os.path.join(path_name, 'Train')
test_copy_path = os.path.join(path_name, 'Test')


for train_dir in train_dirs:
    shutil.move(os.path.join(path_name, train_dir), train_copy_path)

for test_dir in test_dirs:
    shutil.move(os.path.join(path_name, test_dir), test_copy_path)
