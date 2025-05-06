#include <GL/glut.h>
#include <GL/freeglut_ext.h> 
#include <stdlib.h>  // for exit()
#include "viewer.h"

static SimpleGraph* g_graph = NULL;
static float zoom_factor = 10.0f;
static float pan_x = 0.0f, pan_y = 0.0f;
static int is_dragging = 0;
static int last_x, last_y;

void display_callback() {
    if (!g_graph) {
        printf("Graph is NULL!\n");
        return;
    }

    // Projection 설정
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    float range = 10.0f * zoom_factor;
    gluOrtho2D(-range + pan_x, range + pan_x, -range + pan_y, range + pan_y);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glClear(GL_COLOR_BUFFER_BIT);
    glPointSize(5.0f);

    glColor3f(1.0f, 1.0f, 1.0f);
    // glBegin(GL_POINTS);
    for (int i = 0; i < g_graph->num_vertices; ++i) {
        VertexSE2* v = &g_graph->vertices[i];
        // glVertex2f(v->x, v->y);
        glLineWidth(4.0f);
        glBegin(GL_LINES);
        float len = 0.5f;  // 방향선 길이
        float dx = len * cosf(v->theta);
        float dy = len * sinf(v->theta);
        glVertex2f(v->x, v->y);         // 시작점
        glVertex2f(v->x + dx, v->y + dy);  // 끝점

        glEnd();
    }
    glEnd();

    glLineWidth(1.0f);
    glBegin(GL_LINES);
    for (int i = 0; i < g_graph->num_edges; ++i) {
        EdgeSE2* e = &g_graph->edges[i];
        VertexSE2* v1 = &g_graph->vertices[e->from];
        VertexSE2* v2 = &g_graph->vertices[e->to];

        // 정보 행렬의 trace 계산
        float trace = e->information[0][0] + e->information[1][1] + e->information[2][2];
        float intensity = fminf(trace / 100.0f, 1.0f);  // 스케일 조절 (0 ~ 1)

        // 예: 낮은 정보는 빨강, 높은 정보는 초록
        glColor3f(1.0f - intensity, intensity, 0.0f);

        glVertex2f(v1->x, v1->y);
        glVertex2f(v2->x, v2->y);
    }
    glEnd();

    glutSwapBuffers();
}


void reshape_callback(int width, int height) {
    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    float range = 10.0f * zoom_factor;
    gluOrtho2D(-range + pan_x, range + pan_x, -range + pan_y, range + pan_y);}

// 🖱 마우스 클릭 시 프로그램 종료
void mouse_callback(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            is_dragging = 1;
            last_x = x;
            last_y = y;
        } else {
            is_dragging = 0;
        }
    }
    // 마우스 휠 줌 처리 (wheel up = 3, wheel down = 4 in freeglut)
    if (state == GLUT_DOWN) {
        if (button == 3) zoom_factor *= 0.9f;  // zoom in
        if (button == 4) zoom_factor *= 1.1f;  // zoom out\
        
        if (zoom_factor < 0.1f) zoom_factor = 0.1f;
        if (zoom_factor > 100.0f) zoom_factor = 100.0f;
        glutPostRedisplay();
    }
}

void motion_callback(int x, int y) {
    if (!is_dragging) return;

    int dx = x - last_x;
    int dy = y - last_y;
    last_x = x;
    last_y = y;

    // 스크린 좌표 -> 월드 좌표 보정
    float scale = 20.0f * zoom_factor / 800.0f;  // assuming 800x600 window
    pan_x -= dx * scale;
    pan_y += dy * scale;

    glutPostRedisplay();
}

void launch_viewer(SimpleGraph* graph) {
    int argc = 1;
    char* argv[] = {(char*)"viewer"};

    g_graph = graph;

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(800, 600);
    glutCreateWindow("SimplePGO Viewer");

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-10, 10, -10, 10);
    
    glutDisplayFunc(display_callback);
    glutMotionFunc(motion_callback); // 드래그 핸들러 등록
    glutMouseFunc(mouse_callback);  // 👈 클릭 핸들러 등록
    glutReshapeFunc(reshape_callback); // 화면 크기 변경 및 zoom 처리

    glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION);

    glutMainLoop();

}
