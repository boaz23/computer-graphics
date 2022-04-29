#pragma once   //maybe should be static class
#include "igl/opengl/glfw/Display.h"
#include "igl/opengl/glfw/renderer.h"
#include "Assignment2.h"
#include "imgui/imgui.h"


	void glfw_mouse_callback(GLFWwindow* window,int button, int action, int mods)
	{	
		if (action == GLFW_PRESS)
		{
			Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
			Assignment2* scn = (Assignment2*)rndr->GetScene();
			double x2, y2;
			
			glfwGetCursorPos(window, &x2, &y2);
			rndr->UpdatePress(x2, y2);
			scn->intersection(Eigen::Vector3f(x2 / 600 - 1, 1 - y2 / 400, 0));
			if (rndr->Picking((int)x2, (int)y2))
			{
				rndr->UpdatePosition(x2, y2);
				scn->UpdatePosition(x2, y2);
				if(button == GLFW_MOUSE_BUTTON_LEFT)
					rndr->Pressed();
			}
			else
			{
				rndr->UnPick(2);
			}

		}
		else
		{
			Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
			Assignment2* scn = (Assignment2*)rndr->GetScene();
			rndr->UnPick(2);
		}
	}
	
	void glfw_scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
	{
		Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
		Assignment2* scn = (Assignment2*)rndr->GetScene();
		scn->ChangeZoomBy((float)yoffset);
	}
	
	void glfw_cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
	{
		Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
		Assignment2* scn = (Assignment2*)rndr->GetScene();

		rndr->UpdatePosition((float)xpos,(float)ypos);
		scn->UpdatePosition((int)xpos, (int)ypos);

		if (rndr->CheckViewport(xpos,ypos, 0))
		{
			if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS)
			{
				//rndr->MouseProccessing(GLFW_MOUSE_BUTTON_LEFT);
				scn->WhenTranslate();
			}
			else if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS)
			{
				//rndr->MouseProccessing(GLFW_MOUSE_BUTTON_LEFT);
				scn->TransformObject();
			}
		}
	}

	void glfw_window_size_callback(GLFWwindow* window, int width, int height)
	{
		Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);

        rndr->resize(window,width,height);
		
	}
	
	void glfw_key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
	{
		Renderer* rndr = (Renderer*)glfwGetWindowUserPointer(window);
		Assignment2* scn = (Assignment2*)rndr->GetScene();
		//rndr->FreeShapes(2);
		if (action == GLFW_PRESS || action == GLFW_REPEAT)
		{
			switch (key)
			{
			case GLFW_KEY_ESCAPE:
				glfwSetWindowShouldClose(window, GLFW_TRUE);
				break;
				
			case GLFW_KEY_SPACE:
				if (scn->IsActive())
					scn->Deactivate();
				else
					scn->Activate();
				break;

			case GLFW_KEY_UP:
				scn->RotateScene(1, 0);
				break;
			case GLFW_KEY_DOWN:
				scn->RotateScene(-1, 0);
				break;
			case GLFW_KEY_LEFT:
				scn->RotateScene(0, -1);
				break;
			case GLFW_KEY_RIGHT:
				scn->RotateScene(0, 1);
				break;
			case GLFW_KEY_U:
				rndr->MoveCamera(0, scn->yTranslate, 0.25f);
				break;
			case GLFW_KEY_D:
				rndr->MoveCamera(0, scn->yTranslate, -0.25f);
				break;
			case GLFW_KEY_L:
				rndr->MoveCamera(0, scn->xTranslate, -0.25f);
				break;
			
			case GLFW_KEY_R:
				rndr->MoveCamera(0, scn->xTranslate, 0.25f);
				break;
			
			case GLFW_KEY_B:
				rndr->MoveCamera(0, scn->zTranslate, 0.5f);
				break;
			case GLFW_KEY_F:
				rndr->MoveCamera(0, scn->zTranslate, -0.5f);
				break;
			default:
				break;

			}
		}
	}


void Init(Display& display, igl::opengl::glfw::imgui::ImGuiMenu *menu)
{
    display.AddKeyCallBack(glfw_key_callback);
    display.AddMouseCallBacks(glfw_mouse_callback, glfw_scroll_callback, glfw_cursor_position_callback);
    display.AddResizeCallBack(glfw_window_size_callback);
    menu->init(&display);
}
