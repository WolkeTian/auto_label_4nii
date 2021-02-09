# auto_label_4nii
Auto get anatomy label for nifti file, mainly help to identify ICA map.
Support 3D Nifti & 4D Nifti Data
使用：
Results = NII_getlabel('ICA.Maps.nii'); % 阈值默认2
Results = NII_getlabel('ICA.Maps.nii'， threshold); % 手动规定阈值

结果：
报告参考模板（atlas/Networks)，得到每个Region/Network阈上体素的voxel size（数量），百分比（该region阈上体素数量比所有阈上体素数量），和voxel mass（阈上体素的值的和，不止考虑数量也考虑体素值本身大小），以及阈上voxel size最大的区域（MaxRegion)。

注意：
报告的voxel size和voxel mass全部基于atlas分辨率(1mm*1mm*1mm), 除了基于yeo 7/17 Networks报告的的结果（基于输入的nii图像原始分辨率）
