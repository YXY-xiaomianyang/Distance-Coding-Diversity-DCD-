# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 16:24:30 2025

@author: 26744
"""

import numpy as np
from osgeo import gdal

def read_img(filename):
    dataset = gdal.Open(filename)  
    im_width = dataset.RasterXSize  
    im_height = dataset.RasterYSize  
    im_geotrans = dataset.GetGeoTransform()  
    im_proj = dataset.GetProjection()  
    im_data = dataset.ReadAsArray(0, 0, im_width, im_height).astype(np.float64)  
    del dataset  # 关闭对象，文件dataset
    return im_proj, im_geotrans, im_data, im_height, im_width

def write_img(filename, im_proj, im_geotrans, im_data):
    if 'int8' in im_data.dtype.name:
        datatype = gdal.GDT_Byte
    elif 'int16' in im_data.dtype.name:
        datatype = gdal.GDT_UInt16
    else:
        datatype = gdal.GDT_Float32

    if len(im_data.shape) == 3:
        im_bands, im_height, im_width = im_data.shape
    else:
        im_bands, (im_height, im_width) = 1, im_data.shape
        
    driver = gdal.GetDriverByName("GTiff")  
    dataset = driver.Create(filename, im_width, im_height, im_bands, datatype)

    dataset.SetGeoTransform(im_geotrans)  
    dataset.SetProjection(im_proj)  

    if im_bands == 1:
        dataset.GetRasterBand(1).WriteArray(im_data)  
    else:
        for i in range(im_bands):
            dataset.GetRasterBand(i + 1).WriteArray(im_data[i])
    del dataset


def DCD_window(windw, pattern='H', div_mode='Shannon-Weiner'):
    """
    输入一个窗口，计算其在指定连接方式下的距离编码多样性特征值

    Parameters
    ----------
    windw : ndarray
        输入窗口。当pattern参数为'Circle'时，仅能接收窗口大小为3、5或7的窗口

    pattern : str
        连接方式，可选值如下：
            - 'H' : 按水平(Horizontal)方式连接（默认）
            - 'HS' : 按水平蛇形(Horizontal S-shaped)方式连接
            - 'V' : 按竖直(Vertical)方式连接
            - 'VS' : 按竖直蛇形(Vertical S-shaped)方式连接
            - 'Circle' : 按圆形方式连接

    div_mode : str
        多样性度量方式，可选值如下：
            - 'Shannon-Weiner' : 以香农-维纳指数度量多样性（默认）
            - 'Simpson' : 以辛普森指数度量多样性

    Returns
    ----------
    result : float
        返回该窗口的距离编码多样性特征值
    """
    # 如果输入的pattern参数不符合要求，直接给出错误提示
    if pattern not in ['H', 'HS', 'V', 'VS', 'Circle']:
        raise ValueError("pattern参数应为'H'、'HS'、'V'、'VS'或'Circle'")

    # 如果输入的div_mode参数不符合要求，直接给出错误提示
    if div_mode not in ['Shannon-Weiner', 'Simpson']:
        raise ValueError("div_mode参数应为'Shannon-Weiner'或'Simpson'")

    # 按指定的连接方式展平数组
    match pattern:
        case 'H':
            windw_ravel = windw.ravel()   # 按水平跳跃方式展平数组
        case 'HS':
            out = windw.copy()     # 复制窗口
            out[1::2] = windw[1::2, ::-1]  # 利用切片操作将奇数行的灰度值顺序反转
            windw_ravel = out.ravel()      # 展平数组（此时相当于按水平蛇形方式展平）
        case 'V':
            windw = windw.transpose()   # 转置数组（相当于水平转竖直）
            windw_ravel = windw.ravel()   # 展平数组(此时相当于按竖直跳跃方式展平)
        case 'VS':
            windw = windw.transpose()   # 转置数组（相当于水平转竖直）
            out = windw.copy()
            out[1::2] = windw[1::2, ::-1]   # 奇数行灰度值反转（相当于蛇形索引）
            windw_ravel = out.ravel()   # 展平数组（此时相当于竖直蛇形索引）
        case 'Circle':
            win_s = windw.shape[0]  # 窗口大小
            # 如果输入的窗口大小不符合要求，直接给出错误提示
            if win_s not in [3, 5, 7]:
                raise ValueError("pattern参数为'Circle'时，输入窗口的大小应为3、5或7")

            a = windw.ravel()
            # 按圆形排列方式逐个索引灰度值
            match win_s:
                # 窗口大小为3时
                case 3:
                    windw_ravel = np.array([
                        a[4], a[5], a[2], a[1], a[0], a[3], a[6],
                        a[7], a[8]
                    ])
                # 窗口大小为5时
                case 5:
                    windw_ravel = np.array([
                        a[12], a[13], a[14], a[19], a[18], a[24],
                        a[23], a[17], a[22], a[21], a[16], a[20], a[15],
                        a[11], a[10], a[5], a[6], a[0], a[1], a[7],
                        a[2], a[3], a[8], a[4], a[9]
                    ])
                # 窗口大小为7时
                case 7:
                    windw_ravel = np.array([
                        a[24], a[25], a[26], a[27], a[34], a[33], a[41],
                        a[48], a[40], a[32], a[47], a[39], a[46], a[45],
                        a[38], a[31], a[44], a[37], a[43], a[42], a[36],
                        a[30], a[35], a[29], a[28], a[21], a[22], a[23],
                        a[14], a[7], a[15], a[16], a[8], a[0], a[1], a[9],
                        a[2], a[3], a[10], a[17], a[4], a[11], a[5], a[6],
                        a[12], a[18], a[13], a[19], a[20]
                    ])

    # 使用diversity函数计算展平数组的距离编码多样性特征值
    result = diversity(windw_ravel, div_mode)
    return result


def diversity(arr, div_mode='Shannon-Weiner'):
    """
    窗口按某种连接方式展平为一维数组后，计算该数组的距离编码多样性特征值

    Parameters
    ----------
    arr : ndarray
        输入数组，是窗口按某种连接方式被展平后得到的一维数组

    div_mode : str
        多样性度量方式，可选值如下：
            - 'Shannon-Weiner' : 以香农-维纳指数度量多样性（默认）
            - 'Simpson' : 以辛普森指数度量多样性

    Returns
    ----------
    result : float
        输入数组的距离编码多样性特征值
    """

    temp = []  # 创建空列表，用于记录每种灰度值下每种距离编码值出现的频次
    # 第一层循环，遍历所有灰度值
    key = np.unique(arr)  # 找出唯一灰度值
    for i in key:
        indx_tuple = np.where(arr == i)  # 统计灰度值i的索引值
        indx = indx_tuple[0]  # np.where()返回的是tuple，里面包含的是索引值数组
        indx = np.append(indx, indx[0] + arr.size)  # 补充索引值
        diff_indx = np.diff(indx)   # 索引值作差，即距离编码

        # 第二层循环，遍历灰度值i下的所有距离编码
        key_indx = np.unique(diff_indx)  # 找出唯一距离编码值
        for j in key_indx:
            v = diff_indx[diff_indx == j].size    # 统计距离编码j出现的频次
            temp.append(v)           # 记录该频次

    temp = np.array(temp)  # list转为numpy array
    p = np.divide(temp, arr.size)   # 借助频次计算频率

    # 根据div_mode参数选择衡量多样性的指数
    match div_mode:
        case 'Shannon-Weiner':
            result = np.sum(-p*np.log(p))   # 计算香农-维纳指数
        case 'Simpson':
            result = 1 - np.sum(p**2)  # 计算辛普森指数

    return result


def cal_DCD(gray, win_s, div_mode='Shannon-Weiner', patterns=['H', 'HS', 'V', 'VS', 'Circle']):
    """
    计算整幅图像的距离编码多样性特征（可指定要计算哪些模式）

    Parameters
    ----------
    gray : ndarray
        输入图像
    win_s : int
        窗口大小（仅能接收3、5、7）
    div_mode : str
        多样性度量方式（'Shannon-Weiner'或'Simpson'）
    patterns : list
        要计算的模式列表，可选'H', 'HS', 'V', 'VS', 'Circle'

    Returns
    ----------
    DCD : ndarray
        返回整幅图像的距离编码多样性特征
    """
    # 验证输入的patterns是否有效
    valid_patterns = ['H', 'HS', 'V', 'VS', 'Circle']
    for p in patterns:
        if p not in valid_patterns:
            raise ValueError(f"无效的模式 '{p}'，应为 {valid_patterns}")
    
    n1, n2 = gray.shape
    half = int((win_s-1)/2)
    
    # 为每种模式创建存储空间
    results = {}
    if 'H' in patterns:
        results['H'] = np.zeros((n1, n2), dtype=np.float32)
    if 'HS' in patterns:
        results['HS'] = np.zeros((n1, n2), dtype=np.float32)
    if 'V' in patterns:
        results['V'] = np.zeros((n1, n2), dtype=np.float32)
    if 'VS' in patterns:
        results['VS'] = np.zeros((n1, n2), dtype=np.float32)
    if 'Circle' in patterns:
        results['Circle'] = np.zeros((n1, n2), dtype=np.float32)

    for i in range(half, n1-half):
        for j in range(half, n2-half):
            window = gray[i-half:i+half+1, j-half:j+half+1]
            
            for pattern in patterns:
                results[pattern][i, j] = DCD_window(window, pattern, div_mode)
    
    # 如果只选择了一种模式，直接返回该模式的结果
    if len(patterns) == 1:
        return results[patterns[0]]
    # 否则返回所有选择模式的平均值
    else:
        return np.mean(list(results.values()), axis=0)

if __name__ == '__main__':
    
    # 设置输入输出文件路径
    in_file = 'test.tif'
    out_img_file = 'test_DCD.tif'

    proj, geotrans, img = read_img(in_file)   # img的形状为（波段数，行数，列数）
    gray = img[0, :, :]  # 以一个波段为例提取DCD征
    n1 = np.size(gray, 0)  # 图像行数
    n2 = np.size(gray, 1)  # 图像列数
    vmax = np.max(gray)       # 图像最大灰度值
    vmin = np.min(gray)       # 图像最小灰度值
    
    win_s = 5    # 设定窗口大小
    ngl = 16     # 设定灰度级数(number of gray level)
    # 计算OTU图像，这里采用灰度级压缩方式，读者也可选择其它方式
    gray = (gray - vmin)/(vmax - vmin) * ngl # 线性拉伸
    gray = np.floor(gray)   # 向下取整
    gray[gray == ngl] = ngl - 1   # 对最大值单独赋值

    # 用户可以选择要计算的模式（这里示例只计算'H'和'V'模式）, 选择5种模式的均值计算灵敏度高但计算时长相对较长
    selected_patterns = ['H', 'V']  # 用户可以根据自己的实际需要修改这里选择需要的模式
    # 计算DCD
    DCD = cal_DCD(gray, win_s=win_s, div_mode='Shannon-Weiner', patterns=selected_patterns)
    
     # 保存数据
    write_img(out_img_file, proj, geotrans, DCD)