import { defineConfig } from 'vite'

export default defineConfig({
  base: '/phgy_311/',
  server: {
    port: 3000,
    open: true
  },
  build: {
    target: 'es2015',
    outDir: 'dist'
  }
})