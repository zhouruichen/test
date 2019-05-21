# GNSS_SatP.py
#  *按给定的N文件，计算所读取O文件的某历元某卫星的在轨位置
#  *将所计算的卫星位置与下载相应的精密星历做比较，
#  看看广播轨道的精度如何

import math

#### 读取数据
def readData():
    # 读取导航电文文件
    fname_18n = r"E:\大三下\GNSS\bjfs0300.18n"  # 读取数据时改具体路径
    data_18n = []
    f_18n = open(fname_18n)
    n_line = f_18n.readline()
    for n_line in f_18n:
        data_18n.append(n_line)
    f_18n.close()
    # -----------------------------------------------------#

    # 读取精密星历文件
    fname_sp3 = r"E:\大三下\GNSS\igs19862.sp3"
    data_sp3 = []
    f_sp3 = open(fname_sp3)
    s_line = f_sp3.readline()
    # print(s_line, end="")
    for s_line in f_sp3:
        data_sp3.append(s_line)
        #print(s_line,end="")
    f_sp3.close()

    return data_18n,data_sp3
   #-----------------------------------------------------#

#### 将读取的数据转化为可操作的数据
def convertData(data_18n,data_sp3):
    ### 转换卫星导航电文数据
    hour_n = []     # 由于该数据是从0点0分0秒开始计算
    min_n = []      # 分钟数
    sat_num_n = []  # 处理卫星号数

    IODE = []     # 数据/星历发布时间
    Crs = []      #
    deltaN = []   # 平均角速度你的改正值deltaN,单位rad/s
    M0 = []       # toe时刻的平近点角,单位rad

    Cuc = []      # 升交角距u = omega + f 的余弦(rad)
    e = []        # toe时刻的轨道偏心率e
    Cus = []
    sqrtA = []    # 表示长半径A的平方根，单位为m^0.5

    toe = []      # TOE星历的参考时刻
    Cic = []      # 轨道倾角i的余弦
    Omega_Toe = []    # toe时的升交点赤经,rad
    Cis = []      #

    i0 = []        # toe时的轨道倾角，单位rad
    Crc = []       # 卫星至地心距离r(m)的余弦
    omega = []     # toe时的近地点角距，单位rad
    Omega_1 = []   # 升交点赤经的变化率，rad/s

    i_1 = []       # 轨道倾角的变化率，rad/s
    GPS_week = []  # GPS周

    t = []         # 卫星电文发送时刻

    for h in range(7,1430,8):
        # 得到时间和卫星号数
        hour_n.append(int(data_18n[h][12:14]))
        min_n.append(int(data_18n[h][15:17]))
        sat_num_n.append(int(data_18n[h][0:2]))

        # 得到广播轨道—1的数据
        IODE.append(float(data_18n[h+1][3:18])*(10**int(data_18n[h+1][19:22])))
        Crs.append(float(data_18n[h+1][22:37])*(10**int(data_18n[h+1][38:41])))
        deltaN.append(float(data_18n[h+1][41:56])*(10**int(data_18n[h+1][57:60])))
        M0.append(float(data_18n[h+1][60:75])*(10**int(data_18n[h+1][76:79])))

        # 得到广播轨道—2的数据
        Cuc.append(float(data_18n[h + 2][3:18]) * (10 ** int(data_18n[h + 2][19:22])))
        e.append(float(data_18n[h + 2][22:37]) * (10 ** int(data_18n[h + 2][38:41])))
        Cus.append(float(data_18n[h + 2][41:56]) * (10 ** int(data_18n[h + 2][57:60])))
        sqrtA.append(float(data_18n[h + 2][60:75]) * (10 ** int(data_18n[h + 2][76:79])))

        # 得到广播轨道—3的数据
        toe.append(float(data_18n[h + 3][3:18]) * (10 ** int(data_18n[h + 3][19:22])))
        Cic.append(float(data_18n[h + 3][22:37]) * (10 ** int(data_18n[h + 3][38:41])))
        Omega_Toe.append(float(data_18n[h + 3][41:56]) * (10 ** int(data_18n[h + 3][57:60])))
        Cis.append(float(data_18n[h + 3][60:75]) * (10 ** int(data_18n[h + 3][76:79])))

        # 得到广播轨道—4的数据
        i0.append(float(data_18n[h + 4][3:18]) * (10 ** int(data_18n[h + 4][19:22])))
        Crc.append(float(data_18n[h + 4][22:37]) * (10 ** int(data_18n[h + 4][38:41])))
        omega.append(float(data_18n[h + 4][41:56]) * (10 ** int(data_18n[h + 4][57:60])))
        Omega_1.append(float(data_18n[h + 4][60:75]) * (10 ** int(data_18n[h + 4][76:79])))

        # 得到广播轨道—5的数据
        i_1.append(float(data_18n[h + 5][3:18]) * (10 ** int(data_18n[h + 5][19:22])))
        GPS_week.append(float(data_18n[h + 5][41:56]) * (10 ** int(data_18n[h + 5][57:60])))

        # 得到广播轨道—7的数据
        t.append(float(data_18n[h + 7][3:18]) * (10 ** int(data_18n[h + 7][19:22])))
    # ----------------------------------------------------------------------------------------#

    ### 转换GPS精密星历数据.sp3
    hour_s = []     # 精密星历中的小时
    min_s = []      # 精密星历中的分钟

    sat_num_s = []  # 精密星历中的卫星号
    # 精密星历对应的坐标值
    X_s = []
    Y_s = []
    Z_s = []

    for hs in range(21,3092,32):
        hour_s.append(int(data_sp3[hs][14:16]))
        min_s.append(int(data_sp3[hs][17:19]))

        for s in range(1,32):
            sat_num_s.append(int(data_sp3[hs+s][2:4]))
            X_s.append(float(data_sp3[hs+s][5:18]))
            Y_s.append(float(data_sp3[hs+s][19:32]))
            Z_s.append(float(data_sp3[hs+s][33:46]))
    # -----------------------------------------------------#

    #将所有参数归到两个列表里
    data_conv_18n = [hour_n, sat_num_n, IODE, Crs, deltaN, M0, Cuc, e, Cus, sqrtA, toe,
                     Cic, Omega_Toe, Cis, i0, Crc,omega,Omega_1, i_1, GPS_week,  t, min_n]    # 导航电文数据
    data_conv_sp3 = [hour_s, min_s, sat_num_s, X_s, Y_s, Z_s]   # 精密星历数据

    return data_conv_18n, data_conv_sp3
   #---------------------------------------------------------------------------------------#

#### 计算卫星位置
def calSatP(data_conv_18n):
    ## 卫星导航电文数据
    hour_n = data_conv_18n[0]     # 由于该数据是从0点0分0秒开始计算
    min_n = data_conv_18n[21]     # 分钟数
    sat_num_n = data_conv_18n[1]  # 处理卫星号数

    IODE = data_conv_18n[2]  # 数据/星历发布时间
    Crs = data_conv_18n[3]  #
    deltaN = data_conv_18n[4]  # 平均角速度你的改正值deltaN,单位rad/s
    M0 = data_conv_18n[5]  # toe时刻的平近点角,单位rad

    Cuc = data_conv_18n[6]  # 升交角距u = omega + f 的余弦(rad)
    e = data_conv_18n[7]  # toe时刻的轨道偏心率e
    Cus = data_conv_18n[8]
    sqrtA = data_conv_18n[9]  # 表示长半径A的平方根，单位为m^0.5

    toe = data_conv_18n[10]  # TOE星历的参考时刻
    Cic = data_conv_18n[11]  # 轨道倾角i的余弦
    Omega_Toe = data_conv_18n[12]  # toe时的升交点赤经,rad
    Cis = data_conv_18n[13]  #

    i0 = data_conv_18n[14]  # toe时的轨道倾角，单位rad
    Crc = data_conv_18n[15]  # 卫星至地心距离r(m)的余弦
    omega = data_conv_18n[16]  # toe时的近地点角距，单位rad
    Omega_1 = data_conv_18n[17]  # 升交点赤经的变化率，rad/s

    i_1 = data_conv_18n[18]  # 轨道倾角的变化率，rad/s
    GAST_week = data_conv_18n[19]  # GPS周

    # t = data_conv_18n[20]  # 卫星电文发送时刻
    # -----------------------------------------------------#
    GM = 3.986005e14  # GM单位为m^3/s^2

    U_sat_num = []
    nn0 = []
    nn1 = []
    Usat_time = int(input("请输入你想得到1月30号什么时刻的卫星坐标（请输入偶数整点)："))
    print("该时刻可以得到坐标的卫星有：",end="")
    if Usat_time%2 == 0 and Usat_time <= 24 :
        for ui in range(len(hour_n)-13):
            if  Usat_time == hour_n[ui]:
                U_sat_num.append(sat_num_n[ui])
                print(sat_num_n[ui],end=" ")
        for ui in range(len(hour_n)-13):
            if  Usat_time == hour_n[ui]:
                nn0.append(int(ui))
        length_nn0 = len(nn0)
        print()

        Usat_num = int(input("请输入你想得到该时刻GPS卫星坐标的编号："))
        for nni in range(nn0[0],nn0[length_nn0-1]):
            if Usat_num == sat_num_n[nni]:
                id = nni
                break
        print("该卫星在导航电文中的序号id=",id)
    #------------------------------------------------------------------------------------#
        n0 = GM ** 0.5 / (sqrtA[id] ** 3)  ## toe时刻平均角速度
        n = n0 + deltaN[id]  ## 观测时刻的平均角速度
        t = toe[id]

        tk = t - toe[id]
        M = M0[id] + n * tk  ## 观测瞬间卫星的平近点角

        ## 开普勒方程计算偏近点角
        E = M
        E0 = 0
        while abs(E0 - E) >= 10e-12:
            E0 = E
            E = M + e[id] * math.sin(E)

        ## 计算真近点角f
        x1 = (1 - e[id] ** 2) ** 0.5 * math.sin(E)
        x2 = (math.cos(E) - e[id])
        f = math.atan2(x1, x2)

        ## 计算升交角距
        u1 = omega[id] + f

        ## 计算摄动改正项sigmaU,sigmaR,sigmaI
        sigmaU = Cuc[id] * math.cos(2 * u1) + Cus[id] * math.sin(2 * u1)
        sigmaR = Crc[id] * math.cos(2 * u1) + Crs[id] * math.sin(2 * u1)
        sigmaI = Cic[id] * math.cos(2 * u1) + Cis[id] * math.sin(2 * u1)

        ## 对u1,r1,i0进行摄动改正
        a = sqrtA[id] * sqrtA[id]
        u = u1 + sigmaU
        r = a * (1 - e[id] * math.cos(E)) + sigmaR
        i = i0[id] + sigmaI + i_1[id] * tk

        ## 计算卫星在轨道面坐标系中的位置
        x = r * math.cos(u)
        y = r * math.sin(u)

        ## 计算观测瞬间升交点的经度L
        omegaE = 7.2921151467e-5  # 地球自转角速度，rad/s
        L = Omega_Toe[id] + (Omega_1[id] - omegaE) * t - Omega_1[id] * toe[id]
        # -----------------------------------------------------------------#

        ## 计算卫星在瞬时地球坐标系中的位置
        X = x * math.cos(L) - y * math.cos(i) * math.sin(L)
        Y = x * math.sin(L) + y * math.cos(i) * math.cos(L)
        Z = y * math.sin(i)

    else:
        print("输入的整点错误，请重新输入偶数整点数！")

    return X,Y,Z,Usat_time,Usat_num
    # ------------------------------------------------#

#### 将计算的结果输出
def output(data_conv_sp3, X, Y, Z,Usat_time,Usat_num):
    ## GPS精密星历数据
    hour_s = data_conv_sp3[0]  # 精密星历中的小时
    min_s = data_conv_sp3[1]   # 精密星历中的分钟

    sat_num_s = data_conv_sp3[2]  # 精密星历中的卫星号
    # 精密星历对应的坐标值
    X_s = data_conv_sp3[3]
    Y_s = data_conv_sp3[4]
    Z_s = data_conv_sp3[5]
    # ------------------------------------------------#

    n = Usat_time*4*31
    for ids in range(96):
        if hour_s[ids] == Usat_time and min_s[ids] == 0:
            for idn in range(n,n+31):
                if sat_num_s[idn] == Usat_num:
                    dX = abs(X - X_s[idn]*1000)
                    dY = abs(Y - Y_s[idn]*1000)
                    dZ = abs(Z - Z_s[idn]*1000)

    print("计算结果的单位皆为:m")
    print("*由卫星导航电文文件计算得到的卫星坐标(X Y Z) = (",X,Y,Z,")")
    print("计算结果与精密星历坐标做差得到误差为(dX dY dZ) =（",dX,dY,dZ,")")

def main():
    data_18n, data_sp3 = readData()
    data_conv_18n, data_conv_sp3 = convertData(data_18n, data_sp3)
    X,Y,Z,Usat_time,Usat_num = calSatP(data_conv_18n)
    output(data_conv_sp3,X,Y,Z,Usat_time,Usat_num)

main()