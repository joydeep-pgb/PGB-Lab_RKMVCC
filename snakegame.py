import pygame
import random

# Initialize Pygame
pygame.init()

# Constants
WIDTH, HEIGHT = 640, 480
CELL_SIZE = 20
FPS = 10

# Colors
BLACK = (0, 0, 0)
WHITE = (255, 255, 255)
RED = (255, 0, 0)
GREEN = (0, 255, 0)

# Create the window
screen = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption('Snake Game')

# Snake attributes
snake = [(200, 200), (210, 200), (220, 200)]
snake_direction = 'RIGHT'

# Food attributes
food_pos = (random.randint(0, (WIDTH - CELL_SIZE) // CELL_SIZE) * CELL_SIZE,
            random.randint(0, (HEIGHT - CELL_SIZE) // CELL_SIZE) * CELL_SIZE)

# Font for the "Game Over" message
game_over_font = pygame.font.Font(None, 36)

# Clock for controlling the game's FPS
clock = pygame.time.Clock()

def draw_snake(snake):
    for pos in snake:
        pygame.draw.rect(screen, GREEN, (pos[0], pos[1], CELL_SIZE, CELL_SIZE))

def draw_food(food_pos):
    pygame.draw.rect(screen, RED, (food_pos[0], food_pos[1], CELL_SIZE, CELL_SIZE))

def move_snake(snake, direction):
    head = list(snake[0])
    if direction == 'UP':
        head[1] -= CELL_SIZE
    elif direction == 'DOWN':
        head[1] += CELL_SIZE
    elif direction == 'LEFT':
        head[0] -= CELL_SIZE
    elif direction == 'RIGHT':
        head[0] += CELL_SIZE
    return head

def collision_with_wall(snake_head):
    return snake_head[0] >= WIDTH or snake_head[0] < 0 or snake_head[1] >= HEIGHT or snake_head[1] < 0

def collision_with_self(snake_head, snake):
    return snake_head in snake[1:]

def show_game_over():
    text = game_over_font.render("Game Over", True, WHITE)
    text_rect = text.get_rect(center=(WIDTH // 2, HEIGHT // 2))
    screen.blit(text, text_rect)
    pygame.display.update()

def main():
    global snake_direction
    global food_pos
    global snake

    game_running = True

    while game_running:
        screen.fill(BLACK)

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                game_running = False
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_UP and snake_direction != 'DOWN':
                    snake_direction = 'UP'
                elif event.key == pygame.K_DOWN and snake_direction != 'UP':
                    snake_direction = 'DOWN'
                elif event.key == pygame.K_LEFT and snake_direction != 'RIGHT':
                    snake_direction = 'LEFT'
                elif event.key == pygame.K_RIGHT and snake_direction != 'LEFT':
                    snake_direction = 'RIGHT'

        if game_running:
            snake_head = move_snake(snake, snake_direction)
            snake.insert(0, tuple(snake_head))

            if collision_with_wall(snake_head) or collision_with_self(snake_head, snake):
                show_game_over()
                pygame.time.wait(2000)  # Display "Game Over" for 2 seconds
                game_running = False

            if snake_head[0] == food_pos[0] and snake_head[1] == food_pos[1]:
                food_pos = (random.randint(0, (WIDTH - CELL_SIZE) // CELL_SIZE) * CELL_SIZE,
                            random.randint(0, (HEIGHT - CELL_SIZE) // CELL_SIZE) * CELL_SIZE)
            else:
                snake.pop()

            draw_food(food_pos)
            draw_snake(snake)

            pygame.display.update()
            clock.tick(FPS)

    pygame.quit()

if __name__ == "__main__":
    main()
