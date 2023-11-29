from selenium import webdriver
from selenium.webdriver.common.by import By
import time

class ElateAutomation:
    def __init__(self, browser_name="chrome", filename="elastic_tensor.dat"):
        self.browser_name = browser_name
        self.filename = filename
    
    def read_matrix_from_file(self):
        with open(self.filename, 'r') as file:
            matrix = []
            for line in file:
                if not line.strip().startswith("#"):
                    values = line.strip().split()
                    matrix.append(values)
        return matrix

    def start_webdriver(self):
        if self.browser_name == "chrome":
            options = webdriver.ChromeOptions()
            options.add_experimental_option("detach", True)
            return webdriver.Chrome(options=options)
        elif self.browser_name == "firefox":
            return webdriver.Firefox()
        elif self.browser_name == "edge":
            return webdriver.Edge()
        elif self.browser_name == "safari":
            return webdriver.Safari()
        else:
            raise ValueError(f"Unsupported browser: {self.browser_name}")

    def paste_matrix_to_elate(self, matrix):
        matrix_str = '\n'.join([' '.join(map(str, row)) for row in matrix])
        
        driver = self.start_webdriver()
        time.sleep(2)
        driver.get("https://progs.coudert.name/elate")

        textarea = driver.find_element(By.ID, "matrixArea")
        textarea.send_keys(matrix_str)
        process_button = driver.find_element(By.XPATH, "//input[@value='Process this matrix']")
        process_button.click()

    def run(self):
        matrix_from_file = self.read_matrix_from_file()
        self.paste_matrix_to_elate(matrix_from_file)
        
        try:
            while True:
                pass
        except KeyboardInterrupt:
            pass

if __name__ == "__main__":
    print("Choose a browser; when done press Ctrl+C: (chrome, firefox, edge, safari)")
    browser_name = input().lower()
    elate = ElateAutomation(browser_name)
    elate.run()

