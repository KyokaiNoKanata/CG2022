#!/usr/bin/env python
# -*- coding:utf-8 -*-
import math


# 本文件只允许依赖math库


def draw_line(p_list, algorithm):
    """绘制线段

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'，此处的'Naive'仅作为示例，测试时不会出现
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    result = []
    if algorithm == 'Naive':
        if x0 == x1:
            for y in range(y0, y1 + 1):
                result.append((x0, y))
        else:
            if x0 > x1:
                x0, y0, x1, y1 = x1, y1, x0, y0
            k = (y1 - y0) / (x1 - x0)
            for x in range(x0, x1 + 1):
                result.append((x, int(y0 + k * (x - x0))))
    elif algorithm == 'DDA':
        if x0 == x1:
            if y0 > y1:
                y0, y1 = y1, y0
            for y in range(y0, y1):
                result.append((int(x0), int(y)))
        else:
            m = (y1 - y0) / (x1 - x0)
            if abs(m) <= 1:
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                y = y0
                for x in range(x0, x1):
                    result.append((int(x), int(y)))
                    y += m
            else:
                if y0 > y1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                x = x0
                for y in range(y0, y1):
                    result.append((int(x), int(y)))
                    x += 1 / m
    elif algorithm == 'Bresenham':
        if x0 == x1:
            if y0 > y1:
                y0, y1 = y1, y0
            for y in range(y0, y1):
                result.append((int(x0), int(y)))
        elif y0 == y1:
            if x0 > x1:
                x0, x1 = x1, x0
            for x in range(x0, x1):
                result.append((int(x), int(y0)))
        else:
            m = (y1 - y0) / (x1 - x0)
            inc = 1 if m > 0 else -1
            if abs(m) <= 1:
                if x0 > x1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                delta_x, delta_y = x1 - x0, y1 - y0
                p = 2 * delta_y - delta_x
                y = y0
                for x in range(x0, x1):
                    result.append((int(x), int(y)))
                    if p * inc >= 0:
                        p = p + 2 * delta_y - 2 * inc * delta_x
                        y = y + inc
                    else:
                        p = p + 2 * delta_y
            else:
                if y0 > y1:
                    x0, y0, x1, y1 = x1, y1, x0, y0
                x0, y0, x1, y1 = y0, x0, y1, x1
                delta_x, delta_y = x1 - x0, y1 - y0
                p = 2 * delta_y - delta_x
                y = y0
                for x in range(x0, x1):
                    result.append((int(y), int(x)))
                    if p * inc >= 0:
                        p = p + 2 * delta_y - 2 * inc * delta_x
                        y = y + inc
                    else:
                        p = p + 2 * delta_y
    return result


def draw_polygon(p_list, algorithm):
    """绘制多边形

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 多边形的顶点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'DDA'和'Bresenham'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    for i in range(len(p_list)):
        line = draw_line([p_list[i - 1], p_list[i]], algorithm)
        result += line
    return result


def draw_ellipse(p_list):
    """绘制椭圆（采用中点圆生成算法）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 椭圆的矩形包围框左上角和右下角顶点坐标
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    x0, y0 = p_list[0]
    x1, y1 = p_list[1]
    result = []
    rx = abs(x1 - x0) / 2
    ry = abs(y1 - y0) / 2
    xc = (x0 + x1) / 2
    yc = (y0 + y1) / 2
    p1 = pow(ry, 2) - pow(rx, 2) * ry + pow(rx, 2) / 4
    x, y = 0, ry
    result.append((int(x), int(y)))
    while pow(ry, 2) * x < pow(rx, 2) * y:
        if p1 < 0:
            x, y = x + 1, y
            p1 = p1 + 2 * pow(ry, 2) * x + pow(ry, 2)
        else:
            x, y = x + 1, y - 1
            p1 = p1 + 2 * pow(ry, 2) * x - 2 * pow(rx, 2) * y + pow(ry, 2)
        result.append((int(x), int(y)))
        result.append((int(-x), int(y)))
        result.append((int(x), int(-y)))
        result.append((int(-x), int(-y)))
    p2 = pow(ry, 2) * pow(x + 1 / 2, 2) + pow(rx, 2) * pow(y - 1, 2) - pow(rx, 2) * pow(ry, 2)
    while y >= 0:
        if p2 > 0:
            x, y = x, y - 1
            p2 = p2 - 2 * pow(rx, 2) * y + pow(rx, 2)
        else:
            x, y = x + 1, y - 1
            p2 = p2 + 2 * pow(ry, 2) * x - 2 * pow(rx, 2) * y + pow(rx, 2)
        result.append((int(x), int(y)))
        result.append((int(-x), int(y)))
        result.append((int(x), int(-y)))
        result.append((int(-x), int(-y)))
    return translate(result, xc, yc)


def draw_curve(p_list, algorithm):
    """绘制曲线

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 曲线的控制点坐标列表
    :param algorithm: (string) 绘制使用的算法，包括'Bezier'和'B-spline'（三次均匀B样条曲线，曲线不必经过首末控制点）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 绘制结果的像素点坐标列表
    """
    result = []
    if algorithm == 'Bezier':
        ...
    elif algorithm == 'B-spline':
        ...
    return result


def translate(p_list, dx, dy):
    """平移变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param dx: (int) 水平方向平移量
    :param dy: (int) 垂直方向平移量
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    return [[int(p[0] + dx), int(p[1] + dy)] for p in p_list]


def rotate(p_list, x, y, r):
    """旋转变换（除椭圆外）

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 旋转中心x坐标
    :param y: (int) 旋转中心y坐标
    :param r: (int) 顺时针旋转角度（°）
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    r = math.radians(r)
    return [[int(x + (x0 - x) * math.cos(r) - (y0 - y) * math.sin(r)),
             int(y + (x0 - x) * math.sin(r) + (y0 - y) * math.cos(r))] for [x0, y0] in p_list]


def scale(p_list, x, y, s):
    """缩放变换

    :param p_list: (list of list of int: [[x0, y0], [x1, y1], [x2, y2], ...]) 图元参数
    :param x: (int) 缩放中心x坐标
    :param y: (int) 缩放中心y坐标
    :param s: (float) 缩放倍数
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1], [x_2, y_2], ...]) 变换后的图元参数
    """
    return [[int(x0 * s + x * (1 - s)), int(y0 * s + y * (1 - s))] for [x0, y0] in p_list]


def area_code(x_min, y_min, x_max, y_max, x, y):
    res = 0b0
    if x < x_min:
        res |= 0b0001
    if x > x_max:
        res |= 0b0010
    if y < y_min:
        res |= 0b0100
    if y > y_max:
        res |= 0b1000
    return res


def clip(p_list, x_min, y_min, x_max, y_max, algorithm):
    """线段裁剪

    :param p_list: (list of list of int: [[x0, y0], [x1, y1]]) 线段的起点和终点坐标
    :param x_min: 裁剪窗口左上角x坐标
    :param y_min: 裁剪窗口左上角y坐标
    :param x_max: 裁剪窗口右下角x坐标
    :param y_max: 裁剪窗口右下角y坐标
    :param algorithm: (string) 使用的裁剪算法，包括'Cohen-Sutherland'和'Liang-Barsky'
    :return: (list of list of int: [[x_0, y_0], [x_1, y_1]]) 裁剪后线段的起点和终点坐标
    """
    (x0, y0), (x1, y1) = p_list[0], p_list[1]
    if x_min > x_max:
        x_min, x_max = x_max, x_min
    if y_min > y_max:
        y_min, y_max = y_max, y_min
    if algorithm == 'Cohen-Sutherland':
        code0 = area_code(x_min, y_min, x_max, y_max, x0, y0)
        code1 = area_code(x_min, y_min, x_max, y_max, x1, y1)
        if code0 == 0 and code1 == 0:
            return p_list
        elif code0 & code1 != 0:
            return []
        else:
            if code0 == 0b0000:
                x0, y0, x1, y1 = x1, y1, x0, y0
                code0, code1 = code1, code0
            while code0 != 0b0000:
                if code0 & 0b1000 != 0:
                    x0 = round((y_max - y0) * (x1 - x0) / (y1 - y0) + x0)
                    y0 = y_max
                elif code0 & 0b0100 != 0:
                    x0 = round((y_min - y0) * (x1 - x0) / (y1 - y0) + x0)
                    y0 = y_min
                elif code0 & 0b0010 != 0:
                    y0 = round((x_max - x0) * (y1 - y0) / (x1 - x0) + y0)
                    x0 = x_max
                elif code0 & 0b0001 != 0:
                    y0 = round((x_min - x0) * (y1 - y0) / (x1 - x0) + y0)
                    x0 = x_min
                code0 = 0b0000
                x0, y0, x1, y1 = x1, y1, x0, y0
                code0, code1 = code1, code0
            return [[int(x0), int(y0)], [int(x1), int(y1)]]
    elif algorithm == 'Liang-Barsky':
        dx, dy = x1 - x0, y1 - y0
        p = [-dx, dx, -dy, dy]
        q = [x0 - x_min, x_max - x0, y0 - y_min, y_max - y0]
        u1 = float(0)
        u2 = float(1)
        for i in range(0, 4):
            if p[i] == 0 and q[i] < 0:
                return []
            if p[i] == 0:
                continue
            r = q[i] / p[i]
            if p[i] < 0:
                u1 = max(u1, r)
            else:
                u2 = min(u2, r)
            if u1 > u2:
                return []
        return [[round(x0 + u1 * (x1 - x0)), round(y0 + u1 * (y1 - y0))],
                [round(x0 + u2 * (x1 - x0)), round(y0 + u2 * (y1 - y0))]]
    return p_list
