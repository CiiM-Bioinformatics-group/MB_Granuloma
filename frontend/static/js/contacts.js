

const header = document.querySelector("header");

    window.addEventListener("scroll", () => {
      const y = window.scrollY;

      if (y > 0 && y < 300) {
        header.classList.add("bg1");
      }

      if (y >= 300) {
        header.classList.add("bg2");
      }

      if (y == 0) {
        header.classList.remove("bg1", "bg2");
      }
    });